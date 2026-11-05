"""Compare trilobit_processor.edr with EDR's own "Seis to SEGY" output, trace by trace (nothing is written).

    python tests/verify_edr.py <chunk folder with decoded files, *_segy, Trilobit_SHT, Trilobit_RCV, header_edr.txt> ...
"""
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from trilobit_processor import edr  # noqa: E402

TRACE_US = 60_000_000


def read_sht(path):
    rows = {l.split()[0]: [int(v) for v in l.split()[1:]] for l in path.read_text().splitlines() if l.strip()}
    t = lambda r: datetime(r[0], 1, 1) + timedelta(days=r[1] - 1, hours=r[2], minutes=r[3], seconds=r[4])
    return t(rows["PS"]), t(rows["PE"])


def read_rcv(path):
    rc = next(l.split("\t") for l in path.read_text().splitlines() if l.startswith("RC"))
    return edr.Receiver(station=int(rc[2]), line=int(rc[3]), x=float(rc[6]), y=float(rc[7]),
                        skew_us=None if rc[4].strip().lower() == "auto" else int(rc[4]))


def verify(folder):
    folder = Path(folder)
    g1 = next(p for p in folder.glob("*.g1") if not p.name.startswith("."))
    run = edr.read_decoded(g1)
    ps, pe = read_sht(folder / "Trilobit_SHT")
    rcv = read_rcv(folder / "Trilobit_RCV")
    lines = (folder / "header_edr.txt").read_text(errors="ignore").splitlines()
    n_samples = TRACE_US // run.dt_us
    first_us = int((ps - edr.EPOCH).total_seconds()) * 1_000_000
    n_shots = int((pe - ps).total_seconds() * 1_000_000) // TRACE_US
    print(f"{folder}: decoded {run.start} -> {run.end}, skew {run.skew_us} us, SHT {ps} -> {pe} = {n_shots} shots, "
          f"RCV station {rcv.station} skew {rcv.skew_us}", flush=True)
    ok = True
    for comp in edr.DECODED_EXTENSIONS:
        ext = {"h": "h1"}.get(comp, comp)
        segy = folder / f"{run.stem}.{ext}_segy"
        size = 240 + 4 * n_samples
        n_edr = (segy.stat().st_size - 3600) // size
        data = np.memmap(run.files[comp], dtype=">i4", mode="r", offset=edr.DECODED_HEADER_BYTES, shape=(run.n_samples,))
        with open(segy, "rb") as f:
            head = f.read(3600)
            file_header_ok = head == edr.text_header(lines) + edr.binary_header(n_samples, run.dt_us)
            bad_header_bytes, bad_samples, bad_traces = set(), 0, []
            for i in range(min(n_shots, n_edr)):
                t_us = first_us + i * TRACE_US
                shift, rest = edr.drift(run, t_us, rcv.skew_us)
                mine_h = np.frombuffer(edr.trace_header(run, comp, rcv, i + 1, t_us, n_samples, rest), np.uint8)
                mine_s = edr.trace_samples(run, data, t_us, n_samples, shift)
                theirs = f.read(size)
                th = np.frombuffer(theirs[:240], np.uint8)
                diff = np.flatnonzero(mine_h != th)
                bad_header_bytes |= {int(d) + 1 for d in diff if not 124 <= d < 130}
                if mine_s.tobytes() != theirs[240:]:
                    bad_samples += 1
                    if len(bad_traces) < 5:
                        bad_traces.append(i)
        good = file_header_ok and not bad_header_bytes and bad_samples == 0 and n_edr == n_shots
        ok &= good
        print(f"   {comp}: EDR traces {n_edr}, mine {n_shots}; file header identical {file_header_ok}; "
              f"trace header bytes differing (besides 125-130): {sorted(bad_header_bytes) or 'none'}; "
              f"traces with different samples: {bad_samples} {bad_traces}", flush=True)
    return ok


if __name__ == "__main__":
    results = [verify(p) for p in sys.argv[1:]]
    print("\nALL IDENTICAL" if all(results) else "\nDIFFERENCES FOUND")
