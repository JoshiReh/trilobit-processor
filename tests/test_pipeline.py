"""Tests for the Trilobit processor.

    python tests/test_pipeline.py <work_dir> [<real_data_station_folder>]

1. synthetic: two stations, two EDRTest208 exports each; station 0002 has the 59 s gap
   between its exports (the 2026 problem). Checks the check, the time placement
   and the Galperin rotation.
2. real data (optional, e.g. /Volumes/NilRodes/zJoshi/1069): the first 3 h of the
   processed x/y/z/h files are turned back into raw g1/g2/g3/h1 files with a broken
   binary header, processed with 1 h segments and compared with the existing
   DailyChunk output of the notebooks.
"""
import math
import shutil
import struct
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from trilobit_processor import cli, core, edr  # noqa: E402

NS, DT_US = 30000, 2000
# binary header words as EDRTest208 writes them (see edr.BINARY_HEADER_WORDS), plus garbage in 3235
BROKEN = {3221: -2768, 3223: -2768, 3227: 256, 3229: 1024, 3231: 256, 3235: 256, 3255: 256}
FIXED = {3217: DT_US, 3219: DT_US, 3221: NS, 3223: NS, 3225: 5, 3227: 1, 3229: 1, 3231: 1, 3235: 0, 3255: 1}


def text_binary_header():
    h = bytearray(3600)
    h[:3200] = b"C 1 SYNTHETIC TRILOBIT TEST".ljust(3200, b" ")
    struct.pack_into(">h", h, 3216, DT_US)
    struct.pack_into(">h", h, 3224, 5)          # IEEE float
    for pos, val in BROKEN.items():
        struct.pack_into(">h", h, pos - 1, val)
    return h


def rcv(stem):
    """Trilobit_RCV values of a synthetic station: skew, serial, station number."""
    station = int(stem[:4])
    return -1000 * station, 2000 + station, station


def trace_header(i, t, skew, serial, station):
    h = bytearray(240)
    h[0:4] = struct.pack(">I", i + 1)
    h[24:28] = struct.pack(">i", skew)
    h[90:92] = struct.pack(">H", serial)
    h[92:94] = struct.pack(">H", station)
    h[114:116] = struct.pack(">H", NS)
    h[116:118] = struct.pack(">H", DT_US)
    core.write_trace_start_time(h, t)
    return h


def edr_log(start, stop, serial):
    f = lambda t: f"(Day {t.timetuple().tm_yday}) {t:%a %b %d %H:%M:%S %Y} (UTC 0)"
    return (f"START LOG\nDate Start: {f(start)}\nDate Stop: {f(stop)}\nUnit Serial Number: {serial}\n"
            f"Seis Sample Interval (microseconds):{DT_US}\n"
            f"Ascii Info: Svalbard 2026 GlaciGas, Fridtjovhamna, Methane Flairs.\nSTOP LOG\nDate Stop: {f(stop + timedelta(seconds=1))}\n")


def write_run(folder, stem, first_trace, n_traces, log_start, log_stop, seed, skew=None, rcv_of=None):
    """rcv_of: station whose Trilobit_RCV values go into the headers (default: from the stem)."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    data = {comp: rng.standard_normal((n_traces, NS)).astype(np.float32) for comp in core.RAW_EXTENSIONS}
    write_segy_data(folder, stem, first_trace, data, skew, rcv_of)
    (folder / (stem + ".txt")).write_text(edr_log(log_start, log_stop, rcv(rcv_of or stem)[1]))
    return data


def write_raw(folder, name, skew):
    """Minimal EDRTest208 .raw file: only the skew at byte offset 214 matters for the check."""
    raw = bytearray(1024)
    raw[214:218] = struct.pack(">i", skew)
    (folder / name).write_bytes(raw)


def write_segy_data(folder, stem, first_trace, data, skew=None, rcv_of=None):
    station_skew, serial, station = rcv(rcv_of or stem)
    skew = station_skew if skew is None else skew
    for comp, ext in core.RAW_EXTENSIONS.items():
        with open(folder / (stem + ext), "wb") as f:
            f.write(text_binary_header())
            for i in range(len(data[comp])):
                f.write(trace_header(i, first_trace + timedelta(seconds=60 * i), skew, serial, station))
                f.write(data[comp][i].astype(">f4").tobytes())


def read_segy(path):
    raw = path.read_bytes()
    ns = struct.unpack(">H", raw[3600 + 114:3600 + 116])[0]
    size = 240 + 4 * ns
    n = (len(raw) - 3600) // size
    traces = np.frombuffer(raw[3600:3600 + n * size], dtype=np.uint8).reshape(n, size)
    data = traces[:, 240:].copy().view(">f4").astype(np.float32)
    times = core.trace_start_times_us(np.frombuffer(traces.tobytes(), dtype=[("header", "u1", 240), ("data", "u1", size - 240)]))
    return raw[:3600], traces[:, :240], data, times


def test_synthetic(work):
    root = work / "synthetic"
    shutil.rmtree(root, ignore_errors=True)
    t0 = datetime(2026, 7, 27, 16, 0)
    a1 = write_run(root / "0001/1stChunk", "0001_1stChunk_seq1", t0, 20, t0, t0 + timedelta(minutes=20), 1)
    a2 = write_run(root / "0001/2ndChunk", "0001_2ndChunk_seq1", t0 + timedelta(minutes=20), 20,
                   t0 + timedelta(minutes=20), t0 + timedelta(minutes=40), 2)
    b1 = write_run(root / "0002/1stChunk", "0002_1stChunk_seq1", t0, 20, t0, t0 + timedelta(minutes=20), 3)
    late = t0 + timedelta(minutes=20, seconds=59)          # EDRTest208 export starts 59 s late, one trace less
    b2 = write_run(root / "0002/2ndChunk", "0002_2ndChunk_seq1", late, 19,
                   t0 + timedelta(minutes=20), t0 + timedelta(minutes=40), 4)
    (root / "0002/2ndChunk/._0002_2ndChunk_seq1.g1_segy").write_bytes(b"\0" * 10)   # macOS junk file
    # station 0003: 2nd export starts 5 min early (overlap, same samples) and its raw chunk/log starts 35 min early
    c1 = write_run(root / "0003/1stChunk", "0003_1stChunk_seq1", t0, 20, t0, t0 + timedelta(minutes=20), 5)
    c2 = write_run(root / "0003/2ndChunk", "0003_2ndChunk_seq1", t0 + timedelta(minutes=15), 25,
                   t0 - timedelta(minutes=20), t0 + timedelta(minutes=40), 6)
    for comp in core.RAW_EXTENSIONS:                                       # overlapping minutes = same data
        c2[comp][:5] = c1[comp][15:]
    write_segy_data(root / "0003/2ndChunk", "0003_2ndChunk_seq1", t0 + timedelta(minutes=15), c2)

    write_raw(root / "0001/2ndChunk", "0001_2ndChunk.raw", rcv("0001")[0])          # matching .raw: no problem
    write_raw(root / "0003", "full0003.raw", rcv("0003")[0])

    scans = {s.station: s for s in core.scan_all(root)}
    assert scans["0001"].ok, core.check_report(scans.values())
    assert scans["0003"].ok, core.check_report(scans.values())
    notes3 = " | ".join(w.text for w in scans["0003"].warnings)
    assert "overlap of -300.000 s" in notes3 and "EDRTest208 log start" in notes3, notes3
    assert not scans["0002"].ok
    texts = " | ".join(p.text for p in scans["0002"].problems)
    notes = " | ".join(w.text for w in scans["0002"].warnings)
    assert "gap of +59.000 s" in texts and "EDRTest208 log" not in texts, texts
    assert "EDRTest208 log start" in notes and "EDRTest208 log stop" in notes, notes
    print("synthetic check:\n" + core.check_report(scans.values()))

    assert cli.main([str(root), "--hours", "0.25"]) == 1          # stops because of 0002
    assert not list(root.glob("*/0.25hourChunks/*.segy"))
    assert cli.main([str(root), "--hours", "0.25", "--force"]) == 0

    # station 0001: 3 segments 15/15/10 traces, every trace = rotated input of the same minute
    files = sorted((root / "0001/0.25hourChunks").glob("0001_x_*.segy"))
    assert [f.name for f in files] == ["0001_x_chunk0000_20260727T160000.segy",
                                      "0001_x_chunk0001_20260727T161500.segy",
                                      "0001_x_chunk0002_20260727T163000.segy"], files
    g = {c: np.concatenate([a1[c], a2[c]]) for c in core.RAW_EXTENSIONS}
    x_expected, y_expected, z_expected = core.galperin_to_xyz(g["g1"], g["g2"], g["g3"])
    expected = dict(x=x_expected, y=y_expected, z=z_expected, h=g["h"])
    for comp in core.OUTPUT_COMPONENTS:
        parts = [read_segy(f) for f in sorted((root / "0001/0.25hourChunks").glob(f"0001_{comp}_*.segy"))]
        data = np.concatenate([p[2] for p in parts])
        times = np.concatenate([p[3] for p in parts])
        assert np.array_equal(data, expected[comp]), comp
        assert np.array_equal(np.diff(times), np.full(39, 60_000_000)), comp
        hdr = parts[0][0]
        found = {pos: struct.unpack(">h", hdr[pos - 1:pos + 1])[0] for pos in FIXED}
        assert found == FIXED, f"{comp}: binary header {found}"

    # station 0002: minute 16:20 = 59 s zeros + first 1 s of the late export; later minutes shifted by 59 s
    _, _, h, times = read_segy(root / "0002/0.25hourChunks/0002_h_chunk0001_20260727T161500.segy")
    k = 5                                                    # 16:15 + 5 min = 16:20
    assert core.to_datetime(times[k]) == datetime(2026, 7, 27, 16, 20)
    assert np.all(h[k][:59 * 500] == 0)
    assert np.array_equal(h[k][59 * 500:], b2["h"][0][:500])
    assert np.array_equal(h[k + 1][:59 * 500], b2["h"][0][500:])
    assert np.array_equal(h[k + 1][59 * 500:], b2["h"][1][:500])
    _, _, h0, _ = read_segy(root / "0002/0.25hourChunks/0002_h_chunk0000_20260727T160000.segy")
    assert np.array_equal(h0, b1["h"][:15])
    # station 0003: overlap removed, 16:00-16:40 continuous and equal to the input
    h3 = np.concatenate([read_segy(f)[2] for f in sorted((root / "0003/0.25hourChunks").glob("0003_h_*.segy"))])
    t3 = np.concatenate([read_segy(f)[3] for f in sorted((root / "0003/0.25hourChunks").glob("0003_h_*.segy"))])
    assert np.array_equal(h3, np.concatenate([c1["h"], c2["h"][5:]]))
    assert np.array_equal(np.diff(t3), np.full(39, 60_000_000))

    report = sorted(root.glob("trilobit_report_*.txt"))[-1].read_text(encoding="utf-8")
    assert "59.000 s without data filled with zeros" in report, report
    print("synthetic: OK")


def test_rcv(work):
    """Wrong Trilobit_RCV values (the 2026 case: 0989 3rdChunk exported with the skew of 0676)."""
    root = work / "rcv"
    shutil.rmtree(root, ignore_errors=True)
    t0 = datetime(2026, 7, 27, 16, 0)
    t1, t2 = t0 + timedelta(minutes=5), t0 + timedelta(minutes=10)
    write_run(root / "0011/1stChunk", "0011_1stChunk_seq1", t0, 5, t0, t1, 1)
    write_run(root / "0011/2ndChunk", "0011_2ndChunk_seq1", t1, 5, t1, t2, 2, skew=rcv("0012")[0])  # skew of 0012
    write_run(root / "0012/1stChunk", "0012_1stChunk_seq1", t0, 10, t0, t2, 3)
    write_run(root / "0013/1stChunk", "0013_1stChunk_seq1", t0, 10, t0, t2, 4)
    write_raw(root / "0013/1stChunk", "0013_1stChunk.raw", 4711)
    log = root / "0013/1stChunk/0013_1stChunk_seq1.txt"
    log.write_text(log.read_text(encoding="utf-8").replace(f"Unit Serial Number: {rcv('0013')[1]}", "Unit Serial Number: 9999"))
    write_run(root / "0014/1stChunk", "0014_1stChunk_seq1", t0, 10, t0, t2, 5)

    scans = {s.station: s for s in core.scan_all(root)}
    text = {k: " | ".join(p.text for p in s.problems) for k, s in scans.items()}
    print("rcv check:\n" + core.check_report(scans.values()))
    assert "skew not the same in all exports" in text["0011"] and "also used by 0012" in text["0011"], text["0011"]
    assert "also used by 0011" in text["0012"], text["0012"]
    assert "but +4711 us in the .raw file" in text["0013"] and "but 9999 in the EDRTest208 log" in text["0013"], text["0013"]
    assert scans["0014"].ok, text["0014"]

    # Trilobit_RCV of 0014 used for the 2nd export of 0018: grouped with 0014, but the names say 0018
    write_run(root / "0018/1stChunk", "0018_1stChunk_seq1", t0, 5, t0, t1, 6)
    write_run(root / "0018/2ndChunk", "0018_2ndChunk_seq1", t1, 5, t1, t2, 7, rcv_of="0014")
    scans = {s.station: s for s in core.scan_all(root)}
    text = " | ".join(p.text for p in scans["0014"].problems)
    assert "station number 14, but the folder/file names say 0018" in text, text
    print("rcv: OK")


def test_any_names(work):
    """One unchunked station in a folder with any name, files with any name: dropped directly or via its parent."""
    parent = work / "names"
    shutil.rmtree(parent, ignore_errors=True)
    t0 = datetime(2026, 7, 27, 16, 0)
    folder = parent / "adbasdbad_0989"
    data = write_run(folder, "whatever", t0, 30, t0, t0 + timedelta(minutes=30), 8, rcv_of="0989")
    for root in (folder, parent):
        stations = core.find_stations(root)
        assert [(s, f) for s, f, _ in stations] == [("0989", folder)], stations
        scans = core.scan_all(root)
        assert scans[0].ok, core.check_report(scans)
    assert cli.main([str(folder), "--hours", "0.25"]) == 0
    files = sorted((folder / "0.25hourChunks").glob("0989_h_*.segy"))
    assert [f.name for f in files] == ["0989_h_chunk0000_20260727T160000.segy", "0989_h_chunk0001_20260727T161500.segy"]
    assert np.array_equal(np.concatenate([read_segy(f)[2] for f in files]), data["h"])
    assert not core.find_stations(work / "names" / "adbasdbad_0989" / "0.25hourChunks")    # output is not input
    print("any names: OK")


def write_decoded(folder, stem, start, minutes, serial, skew, sync0, sync1, seed):
    """Synthetic "Decode Raw" output: 512-byte start-log header + int32 samples per channel."""
    folder.mkdir(parents=True, exist_ok=True)
    epoch = lambda t: int((t - edr.EPOCH).total_seconds())
    h = bytearray(512)
    struct.pack_into(">IHBbII", h, 0, 0xFFFFFFFF, 506, 0x80, 0, epoch(start), epoch(start) + minutes * 60 + 1)
    struct.pack_into(">i", h, 60, serial)
    struct.pack_into(">H", h, 68, DT_US)
    struct.pack_into(">IIi", h, 206, epoch(sync0), epoch(sync1), skew)
    rng = np.random.default_rng(seed)
    n = minutes * 60 * 1_000_000 // DT_US + 500
    data = {}
    for comp, ext in edr.DECODED_EXTENSIONS.items():
        data[comp] = rng.integers(-500_000, 500_000, n, dtype=np.int32)
        with open(folder / (stem + ext), "wb") as f:
            f.write(h)
            f.write(data[comp].astype(">i4").tobytes())
    # aux5: accelerometer every second (-41, -3, 1053), with one saturated value per minute; aux1 every 60 s
    seconds = epoch(start) + np.arange(minutes * 60, dtype=np.int64)
    acc = np.tile(np.array([-41, -3, 1053], np.int16), (len(seconds), 1))
    acc[::60, 2] = 2047
    rec5 = np.zeros(len(seconds), dtype=[("t", ">u4"), ("v", ">i2", 3)])
    rec5["t"], rec5["v"] = seconds, acc
    rec1 = np.zeros(minutes, dtype=[("t", ">u4"), ("v", ">i2")])
    rec1["t"], rec1["v"] = seconds[::60] + 10, 9600 + np.arange(minutes)
    for name, rec in (("aux5", rec5), ("aux1", rec1)):
        (folder / f"{stem}.{name}").write_bytes(bytes(h) + rec.tobytes())
    return data


def test_decoded(work):
    """Decoded files -> app == EDRTest208 "Seis to SEGY" (edr.write_segy, checked against EDRTest208) -> app, byte for byte."""
    a, b = work / "decoded_in", work / "segy_in"
    shutil.rmtree(a, ignore_errors=True)
    shutil.rmtree(b, ignore_errors=True)
    t0 = datetime(2026, 7, 7, 15, 59, 10)                   # decoded data do not start on a full minute
    stations = {42: (-463445, datetime(2026, 7, 5, 13, 37, 1), datetime(2026, 8, 4, 16, 36, 2)),
                43: (101113, datetime(2026, 7, 5, 13, 28, 58), datetime(2026, 8, 5, 13, 26, 23))}
    for number, (skew, s0, s1) in stations.items():
        write_decoded(a / f"full{number:04d}", f"full{number:04d}", t0, 75, 2000 + number, skew, s0, s1, number)
        (a / f"full{number:04d}" / f"full{number:04d}.txt").write_text(
            edr_log(t0, t0 + timedelta(minutes=75), 2000 + number))
        run = edr.read_decoded(a / f"full{number:04d}" / f"full{number:04d}.g1")
        seg = b / f"seg{number:04d}"
        seg.mkdir(parents=True)
        first_us = int((t0 - edr.EPOCH).total_seconds()) * 1_000_000
        for comp in edr.DECODED_EXTENSIONS:                   # = EDRTest208 with Trilobit_SHT from the first sample
            edr.write_segy(run, comp, seg / f"full{number:04d}.{'h1' if comp == 'h' else comp}_segy",
                           edr.Receiver(station=number), first_us, 75, core.EDR_TRACE_US, edr.DEFAULT_HEADER_LINES)

    scans = core.scan_all(a)
    print("decoded check:\n" + core.check_report(scans))
    assert [s.station for s in scans] == ["0042", "0043"] and all(s.ok for s in scans), core.check_report(scans)
    assert all(r.kind == "decoded" for s in scans for r in s.runs)

    # unfinished "Decode Raw": shorter than the recording in the header, files of different length
    c = work / "decoded_unfinished"
    shutil.rmtree(c, ignore_errors=True)
    write_decoded(c / "full0044", "full0044", t0, 10, 2044, 1000, *stations[43][1:], 44)
    with open(c / "full0044" / "full0044.h1", "r+b") as f:
        f.truncate(512 + 4 * 1000)
    text = " | ".join(p.text for p in core.scan_all(c)[0].problems)
    assert "Decode Raw not finished?" in text and "files differ in length" in text, text

    for args in (["--hours", "0.5"],
                 ["--hours", "0.25", "--start", "2026-07-07 16:09:10", "--end", "2026-07-07 16:59:10"]):
        assert cli.main([str(a), *args, "--quiet"]) == 0
        assert cli.main([str(b), *args, "--quiet"]) == 0
        folder = core.segment_label(float(args[1]))
        for number in stations:
            fa = sorted((a / f"full{number:04d}" / folder).glob("*.segy"))
            fb = sorted((b / f"seg{number:04d}" / folder).glob("*.segy"))
            assert [f.name for f in fa] == [f.name for f in fb] and fa, (fa, fb)
            for x, y in zip(fa, fb):
                if "--start" not in args:
                    assert x.read_bytes()[3200:] == y.read_bytes()[3200:], f"{x.name} differs"
                    continue
                # with a start time the decoded route = EDRTest208 with Trilobit_SHT from that start, the SEG-Y
                # here starts at the first sample: only the drift rest of the first trace (177-178) differs
                xa, ya = np.frombuffer(x.read_bytes(), np.uint8), np.frombuffer(y.read_bytes(), np.uint8)
                d = np.flatnonzero(xa[3200:] != ya[3200:]) + 3200      # textual header: filled in, checked below
                assert len(xa) == len(ya) and all(d >= 3600) and {(i - 3600) % 120240 for i in d} <= {176, 177}, x.name
            if "--start" in args:
                assert fa[0].name.endswith("_20260707T160910.segy"), fa[0].name
            text = fa[0].read_bytes()[:3200].decode("cp500")          # 0042_h_...: filled-in EDRTest208 template
            for part in ("C01 CLIENT: UNIS", "PROJECT: Svalbard 2026 GlaciGas", "LOCATION: Fridtjovhamna",
                         "HYDROPHONE H", f"NODE NUMBER: {number:04d}", "RL: 60000 MILLISECONDS"):
                assert part in text, (part, text[:800])
            rows = (a / f"full{number:04d}" / f"{number:04d}_aux_per_minute.csv").read_text(encoding="utf-8").splitlines()
            assert rows[0] == "time_utc,acc_x,acc_y,acc_z,tilt_deg,battery_mV", rows[0]
            tilt = math.degrees(math.atan2(math.hypot(41, 3), 1053))
            assert len(rows) == 1 + 76 and rows[2] == f"2026-07-07 16:00:00,-41,-3,1053,{tilt:.3f},9601", rows[:3]
            print(f"decoded {number:04d} {args}: {len(fa)} files identical to the EDRTest208 SEG-Y route")
    print("decoded: OK")


def test_real(work, station_folder):
    station = station_folder.name
    n_traces = 180
    root = work / "real"
    shutil.rmtree(root, ignore_errors=True)
    src = {c: station_folder / "1stChunk" / f"{station}_1stChunk_{c}.segy" for c in core.OUTPUT_COMPONENTS}
    raw_dir = root / station / "1stChunk"
    raw_dir.mkdir(parents=True)
    stem = f"{station}_1stChunk_seq1"

    loaded = {}
    for c, p in src.items():
        traces = core.open_traces(p, NS)[:n_traces]
        with open(p, "rb") as f:
            header = bytearray(f.read(3600))
        loaded[c] = (header, np.array(traces["header"]), np.array(traces["data"]).astype(np.float32))
    x, y, z = (loaded[c][2] for c in "xyz")
    # inverse Galperin rotation (the matrix is orthonormal -> inverse = transpose)
    g = dict(g1=-0.81650 * x + 0.57735 * z,
             g2=0.40825 * x - 0.70711 * y + 0.57735 * z,
             g3=0.40825 * x + 0.70711 * y + 0.57735 * z,
             h=loaded["h"][2])
    for comp, ext in core.RAW_EXTENSIONS.items():
        header, trace_headers, _ = loaded["h" if comp == "h" else "x"]
        header = bytearray(header)
        for pos, val in BROKEN.items():
            struct.pack_into(">h", header, pos - 1, val)
        with open(raw_dir / (stem + ext), "wb") as f:
            f.write(header)
            for i in range(n_traces):
                f.write(trace_headers[i].tobytes())
                f.write(g[comp][i].astype(">f4").tobytes())
    start = core.to_datetime(core.trace_start_times_us(core.open_traces(src["x"], NS)[:1])[0])
    serial = core.receiver_info(loaded["x"][1][0])[1]
    (raw_dir / (stem + ".txt")).write_text(edr_log(start, start + timedelta(minutes=n_traces), serial))

    assert cli.main([str(root), "--hours", "1"]) == 0
    daily = {c: next((station_folder / "DailyChunk").glob(f"{station}_{c}_chunk0000_*.segy")) for c in core.OUTPUT_COMPONENTS}
    for comp in core.OUTPUT_COMPONENTS:
        d_hdr, d_th, d_data, _ = read_segy_head(daily[comp], n_traces)
        files = sorted((root / station / core.segment_label(1)).glob(f"{station}_{comp}_*.segy"))
        assert len(files) == 3, files
        for k, f in enumerate(files):
            hdr, th, data, _ = read_segy(f)
            found = {pos: struct.unpack(">h", hdr[pos - 1:pos + 1])[0] for pos in FIXED}
            assert found == FIXED, f"{comp}: binary header {found}"
            # the notebooks wrote 3235 = 1 (sweep frequency, a wrong byte position); everything else must match
            assert hdr[:3234] + hdr[3236:] == d_hdr[:3234] + d_hdr[3236:], f"{comp}: file header differs"
            th_cmp, d_th_k = th.copy(), d_th[60 * k:60 * (k + 1)].copy()
            th_cmp[:, 0:8] = 0                                     # trace numbers restart per file
            d_th_k[:, 0:8] = 0
            assert np.array_equal(th_cmp, d_th_k), f"{comp}: trace headers differ"
            ref = d_data[60 * k:60 * (k + 1)]
            if comp == "h":
                assert np.array_equal(data, ref), "h differs"
            else:                                                   # x/y/z went through rotate + inverse
                err = np.max(np.abs(data - ref)) / np.max(np.abs(ref))
                assert err < 1e-4, f"{comp}: relative error {err}"
        print(f"real {station} {comp}: headers identical, data "
              f"{'identical' if comp == 'h' else 'equal within rounding of the inverse rotation'}")
    print("real data: OK")


def read_segy_head(path, n_traces):
    traces = core.open_traces(path, NS)[:n_traces]
    with open(path, "rb") as f:
        header = f.read(3600)
    return header, np.array(traces["header"]), np.array(traces["data"]).astype(np.float32), None


if __name__ == "__main__":
    work = Path(sys.argv[1])
    work.mkdir(parents=True, exist_ok=True)
    test_synthetic(work)
    test_rcv(work)
    test_any_names(work)
    test_decoded(work)
    if len(sys.argv) > 2:
        test_real(work, Path(sys.argv[2]))
