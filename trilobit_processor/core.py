"""Core logic of the Trilobit processor.

Same logic as the cryo-seep notebooks 0-2, as one pipeline:

  check    - is the EDR export of every station one continuous time series?
  fix      - repair the broken EDR binary header, rotate Galperin g1/g2/g3 -> x/y/z
  segment  - cut the continuous recording into segments of N hours with 60 s output traces

`fix` and `segment` run together, output trace by output trace, so memory use stays
small and the EDR input folders are only read, never modified. Every input trace is
placed by the time stamp in its own trace header: a gap in the recording stays a gap
(filled with zeros) instead of being closed up.
"""
from __future__ import annotations

import os
import re
import struct
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from . import edr

TEXT_HEADER_BYTES = 3200
BINARY_HEADER_BYTES = 400
TRACE_HEADER_BYTES = 240
FIRST_TRACE_BYTE = TEXT_HEADER_BYTES + BINARY_HEADER_BYTES

RAW_EXTENSIONS = {"g1": ".g1_segy", "g2": ".g2_segy", "g3": ".g3_segy", "h": ".h1_segy"}
OUTPUT_COMPONENTS = ("x", "y", "z", "h")
TOLERANCE_US = 1000  # 1 ms

EPOCH = datetime(1970, 1, 1)


def to_datetime(us):
    return EPOCH + timedelta(microseconds=int(us))


def to_us(dt):
    return (dt - EPOCH) // timedelta(microseconds=1)


def fmt_s(us):
    return f"{us / 1e6:+.3f} s"


# ------------------------------------------------------------------ SEG-Y access

def read_layout(path):
    """Samples per trace and sample interval from the FIRST TRACE header
    (the binary header of raw EDR files is broken)."""
    with open(path, "rb") as f:
        f.seek(FIRST_TRACE_BYTE)
        h = f.read(TRACE_HEADER_BYTES)
    if len(h) != TRACE_HEADER_BYTES:
        raise ValueError(f"{path.name}: no complete trace header")
    n_samples = struct.unpack(">H", h[114:116])[0]
    dt_us = struct.unpack(">H", h[116:118])[0]
    if n_samples == 0 or dt_us == 0:
        raise ValueError(f"{path.name}: invalid samples per trace ({n_samples}) or sample interval ({dt_us} us)")
    return n_samples, dt_us


def open_traces(path, n_samples):
    """Memory-mapped view of all traces: traces[i]['header'] / traces[i]['data']."""
    trace_size = TRACE_HEADER_BYTES + 4 * n_samples
    n_traces = (path.stat().st_size - FIRST_TRACE_BYTE) // trace_size
    if n_traces <= 0:
        raise ValueError(f"{path.name}: no complete traces")
    dtype = np.dtype([("header", "u1", TRACE_HEADER_BYTES), ("data", ">f4", n_samples)])
    return np.memmap(path, dtype=dtype, mode="r", offset=FIRST_TRACE_BYTE, shape=(n_traces,))


def trace_start_times_us(traces):
    """Start time of every trace (microseconds since 1970), from header bytes 157-172."""
    fields = np.ascontiguousarray(traces["header"][:, 156:172]).view(">u2").astype(np.int64)
    year, doy, hour, minute, second, _, millisec, microsec = fields.T
    days = (year - 1970).astype("datetime64[Y]").astype("datetime64[D]").astype(np.int64) + doy - 1
    return ((((days * 24 + hour) * 60 + minute) * 60 + second) * 1_000_000 + millisec * 1000 + microsec)


def write_trace_start_time(trace_header, start):
    trace_header[156:158] = struct.pack(">H", start.year)
    trace_header[158:160] = struct.pack(">H", int(start.strftime("%j")))
    trace_header[160:162] = struct.pack(">H", start.hour)
    trace_header[162:164] = struct.pack(">H", start.minute)
    trace_header[164:166] = struct.pack(">H", start.second)
    trace_header[166:168] = struct.pack(">H", 2)  # GMT/UTC
    trace_header[168:170] = struct.pack(">H", start.microsecond // 1000)
    trace_header[170:172] = struct.pack(">H", start.microsecond % 1000)


def fix_binary_header(header, n_samples, dt_us):
    """Patch the fields EDR leaves as garbage (positions as in 1_Binary_Header_Correction)."""
    header = bytearray(header)
    fixes = {3217: dt_us, 3221: n_samples, 3223: dt_us, 3227: n_samples, 3229: 1, 3231: 1, 3235: 1, 3255: 1}
    for byte_1_based, value in fixes.items():
        header[byte_1_based - 1:byte_1_based + 1] = struct.pack(">H", value)
    return header


def receiver_info(trace_header):
    """Values EDR copies from Trilobit_RCV into every trace header: skew [us], unit serial, station number."""
    h = bytes(trace_header)
    return struct.unpack(">i", h[24:28])[0], struct.unpack(">H", h[90:92])[0], struct.unpack(">H", h[92:94])[0]


def read_raw_skew(folder):
    """Skew stored in an EDR .raw file at byte offset 214 (as get_skew.ps1); None if there is no .raw file."""
    raws = sorted(p for p in folder.glob("*.raw") if not p.name.startswith((".", "old_")))
    if not raws:
        return None
    with open(raws[0], "rb") as f:
        f.seek(214)
        b = f.read(4)
    return struct.unpack(">i", b)[0] if len(b) == 4 else None


def galperin_to_xyz(g1, g2, g3):
    x = -0.81650 * g1 + 0.40825 * (g2 + g3)
    y = 0.70711 * (g3 - g2)
    z = 0.57735 * (g1 + g2 + g3)
    return x, y, z


# ------------------------------------------------------------------ discovery

@dataclass
class Problem:
    station: str
    where: str
    text: str

    def __str__(self):
        return f"{self.station} {self.where}: {self.text}"


@dataclass
class Run:
    """One EDR export (one recording sequence) = 4 raw files + log."""
    station: str
    folder: Path
    stem: str
    files: dict
    log: Path | None
    n_samples: int = 0
    dt_us: int = 0
    starts_us: np.ndarray | None = None
    log_start_us: int | None = None
    log_stop_us: int | None = None
    log_dt_us: int | None = None
    skew_us: int | None = None          # Trilobit_RCV values EDR wrote into the trace headers
    serial: int | None = None
    station_number: int | None = None
    kind: str = "segy"                  # "segy" = EDR SEG-Y export, "decoded" = EDR "Decode Raw" output
    decoded: edr.DecodedRun | None = None

    @property
    def name(self):
        return self.folder.name

    @property
    def trace_us(self):
        return self.n_samples * self.dt_us

    @property
    def start_us(self):
        return int(self.starts_us[0])

    @property
    def end_us(self):
        return int(self.starts_us[-1]) + self.trace_us


@dataclass
class StationScan:
    station: str
    folder: Path
    runs: list = field(default_factory=list)
    problems: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    @property
    def ok(self):
        return not self.problems

    @property
    def usable(self):
        return [r for r in self.runs if r.starts_us is not None]


def is_visible(path):
    return not path.name.startswith(".")  # skips macOS "._*" files on external drives


STATION_NUMBER = re.compile(r"(?<!\d)\d{4}(?!\d)")  # 4-digit number in a folder or file name


def header_station(path):
    """Station number EDR wrote from Trilobit_RCV into the first trace header (0 if unreadable)."""
    with open(path, "rb") as f:
        f.seek(FIRST_TRACE_BYTE)
        h = f.read(TRACE_HEADER_BYTES)
    return struct.unpack(">H", h[92:94])[0] if len(h) == TRACE_HEADER_BYTES else 0


def path_numbers(root, path):
    """4-digit numbers in the names from the input folder down to the file."""
    names = [root.name, *path.relative_to(root).parts]
    return {int(n) for name in names for n in STATION_NUMBER.findall(name)}


def name_station(root, path):
    """Station number from the names: the 4-digit number closest to the file (file name first), 0 if none."""
    names = [root.name, *path.relative_to(root).parts]
    for name in reversed(names):
        found = STATION_NUMBER.findall(name)
        if found:
            return int(found[-1])
    return 0


def find_stations(root):
    """All exports below root, grouped by station: [(station, folder, runs)].

    An export is either EDR SEG-Y (*.g1_segy, .g2_segy, .g3_segy, .h1_segy) or the output of EDR
    "Decode Raw" (*.g1, .g2, .g3, .h1); decoded files are preferred if both exist with the same name.
    Folder and file names do not matter: the station number of SEG-Y exports comes from the trace
    headers (only if it is missing there, from a 4-digit number in the names); decoded files have no
    station number, so it comes from the names (e.g. full0989.g1, 0989/...). The station folder,
    where the output goes: the deepest folder containing all exports of the station, or above it
    (up to root) the deepest one whose name contains the station number (e.g. 1069/ for 1069/1stChunk)."""
    root = Path(root)
    groups = {}
    decoded = set()
    for g1 in sorted(root.rglob("*" + edr.DECODED_EXTENSIONS["g1"])):
        if not is_visible(g1):
            continue
        stem = g1.name[:-len(edr.DECODED_EXTENSIONS["g1"])]
        station = f"{name_station(root, g1):04d}"
        files = {c: g1.with_name(stem + ext) for c, ext in edr.DECODED_EXTENSIONS.items()}
        log = g1.with_name(stem + ".txt")
        groups.setdefault(station, []).append(Run(station=station, folder=g1.parent, stem=stem, files=files,
                                                  log=log if log.exists() else None, kind="decoded"))
        decoded.add((g1.parent, stem))
    for g1 in sorted(root.rglob("*" + RAW_EXTENSIONS["g1"])):
        if not is_visible(g1):
            continue
        stem = g1.name[:-len(RAW_EXTENSIONS["g1"])]
        if (g1.parent, stem) in decoded:      # EDR SEG-Y made from decoded files that are used directly
            continue
        number = header_station(g1) or name_station(root, g1)
        station = f"{number:04d}"
        files = {c: g1.with_name(stem + ext) for c, ext in RAW_EXTENSIONS.items()}
        log = g1.with_name(stem + ".txt")
        groups.setdefault(station, []).append(Run(station=station, folder=g1.parent, stem=stem, files=files,
                                                  log=log if log.exists() else None))
    stations = []
    for station, runs in sorted(groups.items()):
        common = Path(os.path.commonpath([r.folder for r in runs]))
        upward = [common, *(p for p in common.parents if p == root or root in p.parents)]
        named = [p for p in upward if int(station) in {int(n) for n in STATION_NUMBER.findall(p.name)}]
        stations.append((station, named[0] if named else common, runs))
    return stations


def read_edr_log(path):
    """Recording start/stop and sample interval from the EDR log (*_seq1.txt)."""
    text = path.read_text(errors="ignore")

    def date(key):
        m = re.search(rf"{key}: \(Day \d+\) \w+ (\w+ +\d+ [\d:]+ \d+)", text)
        return None if m is None else to_us(datetime.strptime(m.group(1), "%b %d %H:%M:%S %Y"))

    m = re.search(r"Seis Sample Interval \(microseconds\):\s*(\d+)", text)
    return date("Date Start"), date("Date Stop"), None if m is None else int(m.group(1))


def read_log_serial(path):
    m = re.search(r"Unit Serial Number:\s*(\d+)", path.read_text(errors="ignore"))
    return None if m is None else int(m.group(1))


# ------------------------------------------------------------------ step 0: check

EDR_TRACE_US = 60_000_000      # EDR record length for decoded exports (RL 60000 ms in Trilobit_RCV)


def set_decoded_grid(run, anchor_us=None, end_us=None):
    """Trace start times of a decoded export, as EDR "Seis to SEGY" with Trilobit_SHT from anchor_us
    (default: first sample) to end_us (default: last sample): 60 s traces, complete ones inside the data."""
    d = run.decoded
    data_start = d.start_epoch * 1_000_000
    data_end = data_start + d.n_samples * d.dt_us
    anchor = data_start if anchor_us is None else anchor_us
    first = anchor + max(0, -(-(data_start - anchor) // run.trace_us)) * run.trace_us
    last_end = data_end if end_us is None else min(end_us, data_end)
    n = max(0, (last_end - first) // run.trace_us)
    run.starts_us = first + run.trace_us * np.arange(n, dtype=np.int64) if n else None


def _scan_segy(scan, run, where, log):
    """EDR SEG-Y export: trace times, layout, gaps inside the files. False if the export is unusable."""
    empty = [p.name for p in run.files.values() if p.stat().st_size <= FIRST_TRACE_BYTE]
    if empty:
        scan.problems.append(Problem(scan.station, where,
                                     f"no traces, only the file header - EDR export failed? ({', '.join(empty)})"))
        return False
    try:
        layouts = {c: read_layout(p) for c, p in run.files.items()}
        if len(set(layouts.values())) > 1:
            scan.problems.append(Problem(scan.station, where, f"components differ in samples/interval: {layouts}"))
            return False
        n_samples, dt_us = layouts["g1"]
        times = {}
        for c, p in run.files.items():
            traces = open_traces(p, n_samples)
            log(f"      reading the times of {len(traces)} traces: {p.name}")
            times[c] = trace_start_times_us(traces)
    except ValueError as e:              # truncated or broken file: report it, check the other runs
        scan.problems.append(Problem(scan.station, where, str(e)))
        return False
    run.n_samples, run.dt_us = n_samples, dt_us

    if any(len(t) != len(times["g1"]) or np.any(t != times["g1"]) for t in times.values()):
        counts = {c: len(t) for c, t in times.items()}
        scan.problems.append(Problem(scan.station, where,
                                     f"components g1/g2/g3/h1 differ in trace times (traces: {counts})"))
    run.starts_us = times["g1"]

    sizes = {c: (p.stat().st_size - FIRST_TRACE_BYTE) % (TRACE_HEADER_BYTES + 4 * run.n_samples)
             for c, p in run.files.items()}
    for c, rest in sizes.items():
        if rest:
            scan.warnings.append(Problem(scan.station, where, f"{c}: {rest} bytes after the last complete trace ignored"))

    steps = np.diff(run.starts_us)
    bad = np.flatnonzero(np.abs(steps - run.trace_us) > TOLERANCE_US)
    if len(bad):
        examples = ", ".join(f"after {to_datetime(run.starts_us[i])}: {fmt_s(steps[i] - run.trace_us)}"
                             for i in bad[:3])
        scan.problems.append(Problem(scan.station, where,
                                     f"{len(bad)} gap(s)/jump(s) inside the file ({examples})"))
    run.skew_us, run.serial, run.station_number = receiver_info(open_traces(run.files["g1"], run.n_samples)[0]["header"])
    return True


def _scan_decoded(scan, run, where, log):
    """EDR "Decode Raw" output: header (start, skew, serial), length. False if the export is unusable."""
    try:
        d = run.decoded = edr.read_decoded(run.files["g1"])
    except ValueError as e:
        scan.problems.append(Problem(scan.station, where, str(e)))
        return False
    if not d.dt_us:
        scan.problems.append(Problem(scan.station, where, "no sample interval in the decoded header"))
        return False
    # an unfinished "Decode Raw" (still running, aborted, or copied too early): data end before the
    # recording end in the header, files of different length, empty log (EDR writes it at the very end)
    unfinished = []
    data_end = d.start_epoch + d.n_samples * d.dt_us / 1e6
    if data_end < d.stop_epoch - 1:
        unfinished.append(f"data end {to_datetime(int(data_end * 1e6))}, but the recording ends "
                          f"{to_datetime(d.stop_epoch * 1_000_000)} (header)")
    sizes = {p.name: p.stat().st_size for p in run.files.values()}
    if len(set(sizes.values())) > 1:
        unfinished.append(f"files differ in length {sizes}")
    if run.log is not None and run.log.stat().st_size == 0:
        unfinished.append(f"{run.log.name} is empty")
    if unfinished:
        scan.problems.append(Problem(scan.station, where, "Decode Raw not finished? " + "; ".join(unfinished)
                                     + " - only the complete part is used"))
    run.n_samples, run.dt_us = EDR_TRACE_US // d.dt_us, d.dt_us
    run.skew_us, run.serial, run.station_number = d.skew_us, d.serial, int(scan.station)
    log(f"      decoded: {d.start} -> {d.end} UTC, serial {d.serial}, skew {d.skew_us:+d} us "
        f"(GPS syncs {to_datetime(d.sync0 * 1_000_000)} / {to_datetime(d.sync1 * 1_000_000)})")
    if int(scan.station) == 0:
        scan.problems.append(Problem(scan.station, where, "no station number: put the 4-digit station number "
                                                          "into the folder or file name (e.g. full0989)"))
    set_decoded_grid(run)
    if run.starts_us is None:
        scan.problems.append(Problem(scan.station, where, "less than one trace (60 s) of data"))
        return False
    return True


def scan_station(station, station_folder, runs, root, log=None):
    log = log or (lambda text: None)
    scan = StationScan(station=station, folder=station_folder, runs=runs)
    for run in scan.runs:
        where = f"{run.name}/{run.stem}"
        log(f"   {where} ({'decoded' if run.kind == 'decoded' else 'EDR SEG-Y'})")
        missing = [p.name for p in run.files.values() if not p.exists()]
        if missing:
            scan.problems.append(Problem(scan.station, where, f"missing files: {', '.join(missing)}"))
            continue
        if not (_scan_decoded if run.kind == "decoded" else _scan_segy)(scan, run, where, log):
            run.starts_us = None
            continue

        if run.log is None:
            scan.warnings.append(Problem(scan.station, where, "no EDR log (*.txt) - start/stop not compared"))
        else:
            # The log describes the extracted raw chunk, not the SHT range: raw chunks are extracted with
            # margins (EDR needs data before the first sample for the drift correction), so only a note.
            run.log_start_us, run.log_stop_us, run.log_dt_us = read_edr_log(run.log)
            if run.log_start_us is not None and abs(run.start_us - run.log_start_us) > TOLERANCE_US:
                scan.warnings.append(Problem(scan.station, where,
                                             f"first trace {to_datetime(run.start_us)} but EDR log start "
                                             f"{to_datetime(run.log_start_us)} ({fmt_s(run.start_us - run.log_start_us)})"))
            if run.log_stop_us is not None and abs(run.end_us - run.log_stop_us) > TOLERANCE_US:
                scan.warnings.append(Problem(scan.station, where,
                                             f"end of last trace {to_datetime(run.end_us)} but EDR log stop "
                                             f"{to_datetime(run.log_stop_us)} ({fmt_s(run.end_us - run.log_stop_us)})"))
            if run.log_dt_us is not None and run.log_dt_us != run.dt_us:
                scan.problems.append(Problem(scan.station, where,
                                             f"sample interval {run.dt_us} us in trace headers, {run.log_dt_us} us in EDR log"))

        # Trilobit_RCV values (SEG-Y) / header values (decoded): the skew sets the clock drift correction,
        # a wrong one shifts the data
        source = "Trilobit_RCV" if run.kind == "segy" else "decoded header:"
        raw_skew = read_raw_skew(run.folder)
        if raw_skew is None:
            raw_skew = read_raw_skew(station_folder)
        if raw_skew is not None and raw_skew != run.skew_us:
            scan.problems.append(Problem(scan.station, where, f"{source} skew {run.skew_us:+d} us, "
                                                              f"but {raw_skew:+d} us in the .raw file"))
        log_serial = read_log_serial(run.log) if run.log is not None else None
        if log_serial is not None and log_serial != run.serial:
            scan.problems.append(Problem(scan.station, where, f"{source} serial {run.serial}, "
                                                              f"but {log_serial} in the EDR log"))
        # names with 4-digit numbers (e.g. "0989/2ndChunk") must contain this station: catches a
        # Trilobit_RCV copied from another station, whose exports would otherwise be merged with that one
        numbers = path_numbers(root, run.files["g1"])
        if run.kind == "segy" and numbers and run.station_number not in numbers:
            scan.problems.append(Problem(scan.station, where,
                                         f"Trilobit_RCV station number {run.station_number}, but the folder/file "
                                         f"names say {', '.join(f'{n:04d}' for n in sorted(numbers))}"))

    usable = scan.usable
    skews = sorted({r.skew_us for r in usable})
    if len(skews) > 1:
        scan.problems.append(Problem(scan.station, "-", "Trilobit_RCV skew not the same in all exports: "
                                     + ", ".join(f"{r.name} {r.skew_us:+d} us" for r in usable)))
    if len({(r.n_samples, r.dt_us) for r in usable}) > 1:
        scan.problems.append(Problem(scan.station, "-", "runs have different samples per trace / sample interval"))

    usable.sort(key=lambda r: r.start_us)
    scan.runs = usable + [r for r in scan.runs if r.starts_us is None]
    for prev, run in zip(usable[:-1], usable[1:]):
        gap = run.start_us - prev.end_us
        if abs(gap) > TOLERANCE_US:
            # a gap means missing data; an overlap is harmless (every trace is placed by its header time)
            kind, target = ("gap", scan.problems) if gap > 0 else ("overlap", scan.warnings)
            target.append(Problem(scan.station, f"{prev.name} -> {run.name}",
                                  f"{kind} of {fmt_s(gap)} ({to_datetime(prev.end_us)} -> {to_datetime(run.start_us)})"))
    return scan


def scan_all(root, progress=None, log=None):
    log = log or (lambda text: None)
    root = Path(root)
    log(f"Searching for exports in {root} ...")
    stations = find_stations(root)
    log(f"Found {len(stations)} station(s): " + ", ".join(f"{s} ({len(r)} export(s))" for s, _, r in stations))
    scans = []
    for i, (station, folder, runs) in enumerate(stations):
        if progress:
            progress(i / max(len(stations), 1), f"Checking {station} ({i + 1}/{len(stations)}) ...")
        log(f"Checking {station} ({i + 1}/{len(stations)}) in {folder}")
        scans.append(scan_station(station, folder, runs, root, log))
    # the same skew at two stations = Trilobit_RCV copied from another station
    users = {}
    for scan in scans:
        for skew in {r.skew_us for r in scan.usable if r.skew_us}:
            users.setdefault(skew, []).append(scan)
    for skew, same in users.items():
        for scan in same:
            if len(same) > 1:
                others = ", ".join(s.station for s in same if s is not scan)
                scan.problems.append(Problem(scan.station, "-", f"Trilobit_RCV skew {skew:+d} us also used by {others}"))
    return scans


def check_report(scans):
    lines = []
    for scan in scans:
        n_traces = sum(len(r.starts_us) for r in scan.usable)
        status = "OK" if scan.ok else "PROBLEM"
        lines.append(f"[{status}] {scan.station}: {len(scan.runs)} EDR export(s), {n_traces} traces")
        for r in scan.usable:
            kind = " (decoded)" if r.kind == "decoded" else ""
            lines.append(f"      {r.name}/{r.stem}{kind}: {to_datetime(r.start_us)} -> {to_datetime(r.end_us)} UTC, "
                         f"{len(r.starts_us)} traces, skew {r.skew_us:+d} us, serial {r.serial}")
        lines += [f"      PROBLEM {p.where}: {p.text}" for p in scan.problems]
        lines += [f"      note {w.where}: {w.text}" for w in scan.warnings]
    return "\n".join(lines)


# ------------------------------------------------------------------ steps 1+2: fix + segment

@dataclass
class Settings:
    segment_hours: float = 24
    output_trace_seconds: float = 60
    start_utc: datetime | None = None
    end_utc: datetime | None = None


def parse_utc(text):
    """'YYYY-MM-DD HH:MM[:SS]' (UTC) -> datetime; empty -> None."""
    text = (text or "").strip()
    if not text:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass
    raise ValueError(f"'{text}' is not YYYY-MM-DD HH:MM[:SS]")


def segment_label(hours):
    return f"{hours:g}hourChunks"


def output_folder(scan, settings):
    return scan.folder / segment_label(settings.segment_hours)


def existing_outputs(scans, settings):
    return [f for s in scans for f in output_folder(s, settings).glob(f"{s.station}_*_chunk*.segy") if is_visible(f)]


@dataclass
class Plan:
    scan: StationScan
    windows: list
    n_samples: int
    dt_us: int
    out_samples: int

    @property
    def out_trace_us(self):
        return self.out_samples * self.dt_us

    def n_output_traces(self):
        return sum(-(-(b - a) // self.out_trace_us) for a, b in self.windows)


def plan_station(scan, settings):
    # decoded exports: traces like EDR with Trilobit_SHT = start/end of the settings (or of the data)
    for r in scan.runs:
        if r.decoded is not None:
            set_decoded_grid(r, None if settings.start_utc is None else to_us(settings.start_utc),
                             None if settings.end_utc is None else to_us(settings.end_utc))
    runs = scan.usable
    if not runs:
        return None
    n_samples, dt_us = runs[0].n_samples, runs[0].dt_us
    out_samples = int(round(settings.output_trace_seconds * 1e6 / dt_us))
    rec_start = min(r.start_us for r in runs)
    rec_end = max(r.end_us for r in runs)
    start = rec_start if settings.start_utc is None else max(to_us(settings.start_utc), rec_start)
    end = rec_end if settings.end_utc is None else min(to_us(settings.end_utc), rec_end)
    segment_us = int(round(settings.segment_hours * 3600e6))

    starts = np.sort(np.concatenate([r.starts_us for r in runs]))
    windows = []
    t = start
    while t < end:
        t_end = min(t + segment_us, end)
        # keep only windows that contain data (like the notebook: empty chunks are skipped)
        if np.searchsorted(starts, t_end, "left") > np.searchsorted(starts, t - n_samples * dt_us, "right"):
            windows.append((t, t_end))
        else:
            windows.append(None)
        t = t_end
    return Plan(scan=scan, windows=windows, n_samples=n_samples, dt_us=dt_us, out_samples=out_samples)


def processing_summary(plans, settings):
    """Overview before processing (like EDR): sample interval, recording start/end, new start/end, chunks."""
    lines = []
    for p in plans:
        runs = p.scan.usable
        rec_start = min(to_us(r.decoded.start) if r.decoded else r.start_us for r in runs)
        rec_end = max(to_us(r.decoded.end) if r.decoded else r.end_us for r in runs)
        windows = [w for w in p.windows if w is not None]
        lines += [f"Station {p.scan.station}",
                  f"   Sample interval:  {p.dt_us} us ({1e6 / p.dt_us:g} Hz)",
                  f"   Recording start:  {to_datetime(rec_start):%Y-%m-%d %H:%M:%S} UTC",
                  f"   Recording end:    {to_datetime(rec_end):%Y-%m-%d %H:%M:%S} UTC",
                  f"   New start:        {to_datetime(windows[0][0]):%Y-%m-%d %H:%M:%S} UTC" if windows else "   New start: -",
                  f"   New end:          {to_datetime(windows[-1][1]):%Y-%m-%d %H:%M:%S} UTC" if windows else "   New end: -"]
    return "\n".join(lines)


def processing_sentence(plans, settings):
    n = sum(1 for p in plans for w in p.windows if w is not None)
    decoded = any(r.decoded for p in plans for r in p.scan.usable)
    text = (f"The {'decoded ' if decoded else ''}data from EDR will now be chunked into {n} {settings.segment_hours:g}-hour "
            f"chunks per component (x, y, z, h)")
    if decoded:
        text += ", with additional information (tilt, battery, temperature) in a CSV"
    return text + "."


class Cancelled(Exception):
    pass


def header_lines(run):
    """Textual header for decoded exports: header_edr.txt next to the files (as EDR), else EDR's default."""
    for folder in (run.folder, run.folder.parent):
        p = folder / "header_edr.txt"
        if p.exists():
            return p.read_text(errors="ignore").splitlines()
    return edr.DEFAULT_HEADER_LINES


def trace_reader(run):
    """Function (component, index) -> samples of input trace `index` of this export, as in EDR's SEG-Y."""
    if run.decoded is None:
        traces = {c: open_traces(p, run.n_samples) for c, p in run.files.items()}
        return lambda c, i: traces[c][i]["data"]
    d = run.decoded
    data = {c: np.memmap(p, dtype=">i4", mode="r", offset=edr.DECODED_HEADER_BYTES, shape=(d.n_samples,))
            for c, p in d.files.items()}

    def read(c, i):
        t_us = int(run.starts_us[i])
        return edr.trace_samples(d, data[c], t_us, run.n_samples, edr.drift(d, t_us)[0])
    return read


def export_headers(run, src):
    """(3600-byte file header, first trace header) of component src of an export, as in EDR's SEG-Y."""
    if run.decoded is None:
        with open(run.files[src], "rb") as f:
            return f.read(FIRST_TRACE_BYTE), f.read(TRACE_HEADER_BYTES)
    d = run.decoded
    t_us = int(run.starts_us[0])
    file_header = edr.text_header(header_lines(run)) + edr.binary_header(run.n_samples, d.dt_us)
    rcv = edr.Receiver(station=run.station_number)
    return file_header, edr.trace_header(d, src, rcv, 1, t_us, run.n_samples, edr.drift(d, t_us)[1])


def write_aux_csv(scan):
    """<station>_aux_per_minute.csv next to the segment folder (decoded exports only), or None."""
    cols = edr.aux_per_minute([r for r in scan.runs if r.decoded is not None])
    if cols is None:
        return None
    path = scan.folder / f"{scan.station}_aux_per_minute.csv"
    names = list(cols)
    fmt = {"time_utc": lambda v: f"{v:%Y-%m-%d %H:%M:%S}", "tilt_deg": lambda v: f"{v:.3f}"}
    with open(path, "w") as f:
        f.write(",".join(names) + "\n")
        for row in zip(*cols.values()):
            f.write(",".join("" if (not isinstance(v, datetime) and np.isnan(v)) else fmt.get(n, lambda v: f"{v:g}")(v)
                             for n, v in zip(names, row)) + "\n")
    return path


def process_station(plan, settings, progress=None, cancel=None, log=None):
    """Write the segment files of one station. Returns a summary dict."""
    log = log or (lambda text: None)
    scan = plan.scan
    runs = scan.usable
    ns, dt_us, out_ns = plan.n_samples, plan.dt_us, plan.out_samples
    trace_in_us = ns * dt_us
    out_trace_us = plan.out_trace_us

    readers = [trace_reader(r) for r in runs]
    starts = np.concatenate([r.starts_us for r in runs])
    run_of = np.concatenate([np.full(len(r.starts_us), k) for k, r in enumerate(runs)])
    index_of = np.concatenate([np.arange(len(r.starts_us)) for r in runs])
    order = np.argsort(starts, kind="stable")
    starts, run_of, index_of = starts[order], run_of[order], index_of[order]

    # headers: from the earliest export (runs are sorted by start time)
    templates = {}
    for comp, src in (("x", "g1"), ("y", "g1"), ("z", "g1"), ("h", "h")):
        file_header, first_trace_header = export_headers(runs[0], src)
        file_header = fix_binary_header(file_header, ns, dt_us)
        file_header[3220:3222] = struct.pack(">H", out_ns)
        file_header[3222:3224] = struct.pack(">H", out_ns)
        templates[comp] = (bytes(file_header), first_trace_header)

    folder = output_folder(scan, settings)
    folder.mkdir(exist_ok=True)
    summary = dict(station=scan.station, folder=folder, files=0, traces=0, gap_seconds=0.0)
    done_traces, total_traces = 0, plan.n_output_traces()
    used = [w for w in plan.windows if w is not None]
    n_segments = len(used)
    log(f"{scan.station}: {n_segments} segment(s) of {settings.segment_hours:g} h, "
        f"{to_datetime(used[0][0])} -> {to_datetime(used[-1][1])} UTC, output -> {folder}")
    segment_no = 0

    for chunk_index, window in enumerate(plan.windows):
        if window is None:
            continue
        segment_no += 1
        segment_clock = time.monotonic()
        segment_gap = summary["gap_seconds"]
        seg_start, seg_end = window
        n_out = -(-(seg_end - seg_start) // out_trace_us)
        stamp = to_datetime(seg_start)
        paths = {c: folder / f"{scan.station}_{c}_chunk{chunk_index:04d}_{stamp:%Y%m%dT%H%M%S}.segy"
                 for c in OUTPUT_COMPONENTS}
        out = {c: open(p, "wb") for c, p in paths.items()}
        try:
            for c in OUTPUT_COMPONENTS:
                out[c].write(templates[c][0])

            for k in range(n_out):
                if cancel is not None and cancel.is_set():
                    raise Cancelled()
                t0 = seg_start + k * out_trace_us
                n_valid = int(round((min(t0 + out_trace_us, seg_end) - t0) / dt_us))
                buf = {c: np.zeros(out_ns, dtype=np.float32) for c in RAW_EXTENSIONS}
                filled = np.zeros(n_valid, dtype=bool)

                lo = np.searchsorted(starts, t0 - trace_in_us, "right")
                hi = np.searchsorted(starts, t0 + n_valid * dt_us, "left")
                for j in range(lo, hi):
                    offset = int(round((starts[j] - t0) / dt_us))   # position of input sample 0 in this trace
                    a, b = max(0, offset), min(n_valid, offset + ns)
                    if b <= a:
                        continue
                    read = readers[run_of[j]]
                    for c in RAW_EXTENSIONS:
                        buf[c][a:b] = read(c, index_of[j])[a - offset:b - offset]
                    filled[a:b] = True
                summary["gap_seconds"] += (n_valid - int(filled.sum())) * dt_us / 1e6

                x, y, z = galperin_to_xyz(buf["g1"], buf["g2"], buf["g3"])
                for c, samples in zip(OUTPUT_COMPONENTS, (x, y, z, buf["h"])):
                    trace_header = bytearray(templates[c][1])
                    trace_header[0:4] = struct.pack(">I", k + 1)
                    trace_header[4:8] = struct.pack(">I", k + 1)
                    trace_header[114:116] = struct.pack(">H", out_ns)
                    trace_header[116:118] = struct.pack(">H", dt_us)
                    write_trace_start_time(trace_header, to_datetime(t0))
                    out[c].write(trace_header)
                    out[c].write(samples.astype(">f4").tobytes())

                done_traces += 1
                if progress and (done_traces % 10 == 0 or k == n_out - 1):
                    progress(done_traces / total_traces, f"{scan.station}: segment {segment_no}/{n_segments} "
                                                         f"({stamp:%Y-%m-%d %H:%M}), trace {k + 1}/{n_out}")
        except BaseException:
            for f in out.values():
                f.close()
            for p in paths.values():
                p.unlink(missing_ok=True)  # never leave a half-written segment behind
            raise
        for f in out.values():
            f.close()
        summary["files"] += len(paths)
        summary["traces"] += n_out
        gap = summary["gap_seconds"] - segment_gap
        log(f"   segment {segment_no}/{n_segments}: {stamp:%Y-%m-%d %H:%M} -> {to_datetime(seg_end):%Y-%m-%d %H:%M} UTC, "
            f"{n_out} traces x/y/z/h written in {time.monotonic() - segment_clock:.0f} s"
            + (f", {gap:.3f} s without data (zeros)" if gap else ""))
    if any(r.decoded is not None for r in runs):
        log(f"   {scan.station}: writing tilt / battery / aux channels per minute ...")
        summary["aux_csv"] = write_aux_csv(scan)
        if summary["aux_csv"]:
            log(f"   -> {summary['aux_csv']}")
    return summary
