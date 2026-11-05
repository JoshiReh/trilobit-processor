"""EDR "Seis to SEGY" without EDR: decoded files (*.g1/.g2/.g3/.h1 from "Decode Raw") -> SEG-Y.

Reverse-engineered from EDR208 output (station 0676, 2026) and checked against it trace by trace.

Decoded file  = 512-byte header (copy of the start-log packet of the .raw file) + int32 big-endian
                samples, one channel, continuous, first sample at the epoch in header bytes 8-11.
Header bytes  = 8 epoch of the first sample, 12 end of the recording, 60 unit serial, 68 sample interval [us],
                206 last GPS sync before the recording, 210 GPS sync after it ("skew check"),
                214 skew [us] measured at the second sync (same as get_skew.ps1).
Drift         = the skew grows linearly between the two syncs. EDR evaluates it once at the first
                sample of the decoded file (start), truncated to whole us ("MicroSeconds skew before
                first sample"), and adds rate * (t - start) for a trace starting at t, with the rate
                skew / (sync1 - sync0) [us/s] in single precision; the sum is truncated to whole us.
                The trace is the decoded data shifted by trunc(skew(t) / dt) samples (constant within
                the trace), the rest (< 1 sample) goes into trace header bytes 177-178. Samples
                outside the decoded data are zero. (All truncations toward zero.)
"""
from __future__ import annotations

import math
import struct
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

DECODED_HEADER_BYTES = 512
DECODED_EXTENSIONS = {"g1": ".g1", "g2": ".g2", "g3": ".g3", "h": ".h1"}
CHANNEL_CODE = {"h": 11, "g1": 12, "g2": 13, "g3": 14}      # trace header bytes 29-30
EPOCH = datetime(1970, 1, 1)

# EDR's header_edr.txt (40 lines, copied verbatim from the 2026 processing folders, including its
# placeholders); used when there is no header_edr.txt next to the decoded files
DEFAULT_HEADER_LINES = [
    'C01 CLIENT:                                                            UoB/NORCE',
    'C02 PROJECT:     S-TRANET           COUNTRY:                   Nesbyen    NORWAY',
    'C03 CONTRACTOR: SEABED GEOSOLUTIONS    YEAR:                                2026',
    'C04 ============================================================================',
    'C05 INSTRUMENT: SEABED GEOSOLUTIONS TRILOBIT OCEAN BOTTOM NODE (4 COMPONENT )   ',
    'C06 PRODUCT TYPE: RECEIVER GATHER                                   HYDROPHONE H',
    'C07 FORMAT: SEGY IEEE 32BIT SI: 2000 MICROSECONDS  RL:         4000 MILLISECONDS',
    'C08 ============================================================================',
    'C09 NODE NUMBER:    [node_station] EDR SERIAL: [edr]                            ',
    'C10 DATE OF DEPLOYMENT: 2026.04.15 DATE OF RETRIEVAL: 2026.05.20                ',
    'C11 ============================================================================',
    'C12                     SEGY TRACE HEADERS DEFINITION                           ',
    'C13 BYTE 001-004: TRACE NUMBER           BYTE 005-008: TRACE NUMBER             ',
    'C14 BYTE 009-012: NODE STATION NUMBER    BYTE 017-020: SHOT POINT NUMBER        ',
    'C15 BYTE 025-028: CLOCK DRIFT(MICROSEC)  BYTE 029-030: TRACE CODE               ',
    'C16 BYTE 031-032: SOURCE CODE            BYTE 045-048: TIDE VALUE (CM)          ',
    'C17 BYTE 049-052: SOURCE DEPTH (CM)      BYTE 053-056: DEPLOY TIME (EPOCH)      ',
    'C18 BYTE 057-060: RECOVERY TIME (EPOCH)  BYTE 061-064: SOURCE WATER DEPTH (CM)  ',
    'C19 BYTE 065-068: NODE WATER DEPTH (CM)  BYTE 069-070: ELEVATION SCALING FACTOR ',
    'C20 BYTE 071-072: COORD SCALING FACTOR   BYTE 073-076: SOURCE X                 ',
    'C22 BYTE 077-080: SOURCE Y               BYTE 081-084: RECEIVER X               ',
    'C22 BYTE 085-088: RECEIVER Y             BYTE 089-090: COORD UNITS ( 1-> METERS)',
    'C23 BYTE 091-092: EDR NUMBER             BYTE 093-094: TSN NUMBER               ',
    'C24 BYTE 105-106: APPLIED DRIFT (mS)     BYTE 107-108: APPLIED DRIFT (uS)       ',
    'C25 BYTE 115-116: NUMBER OF SAMPLES      BYTE 117-118: SI (MICROSECONDS)        ',
    'C26 BYTE 125-126: ACC Gx (G *10000    )  BYTE 127-128: ACC Gy  (G *10000    )   ',
    'C27 BYTE 129-130: ACC Gz (G *10000    )  BYTE 131-132: HEADING (10TH OF DEG)    ',
    'C28 BYTE 139-140: SHOT SEQUENCE NUMBER   BYTE 157-158: YEAR                     ',
    'C29 BYTE 159-160: DAY OF YEAR            BYTE 161-162: HOUR                     ',
    'C30 BYTE 163-164: MINUTE                 BYTE 165-166: SECOND                   ',
    'C31 BYTE 167-168: TIME BASECODE(2->GMT)  BYTE 169-170: MILLISECONDS             ',
    'C32 BYTE 171-172: MICROSECONDS           BYTE 177-178: CLOCK DRIFT FRACTION     ',
    'C33                     COMMENTS:                                               ',
    'C34 SOURCE CODE: 1/2=good flip/flop, 4=outside positioning spec                 ',
    'C35                7=edited by QC, 8=bad bubble , 9=bad timing                  ',
    'C36 HYDROPHONE SENSITIVITY: 2.81e-4 uBar/bit                                    ',
    'C37 GEOPHONE SENSITIVITY: 8.97e-6 (um/s)/bit                                    ',
    'C38                                                                             ',
    'C39                                                                             ',
    'C40  END OF EBCDIC                                                              ',
]

# EDR's binary header, as written for 60 s records at 2 ms (only these words are non-zero);
# 3217/3219 = sample interval, 3221/3223 = samples per trace with the top bit set
BINARY_HEADER_WORDS = {3207: 1, 3211: 1, 3213: 1, 3225: 5, 3227: 256, 3229: 1024, 3231: 256, 3255: 256}


def _i32(b, o):
    return struct.unpack(">i", b[o:o + 4])[0]


def _u32(b, o):
    return struct.unpack(">I", b[o:o + 4])[0]


@dataclass
class DecodedRun:
    """One "Decode Raw" output: four channel files with the same name."""
    folder: Path
    stem: str
    files: dict
    start_epoch: int          # time of the first sample [s since 1970]
    stop_epoch: int           # end of the recording (or raw chunk) according to the header
    serial: int
    dt_us: int
    sync0: int                # [s since 1970]
    sync1: int
    skew_us: int
    n_samples: int            # per channel
    header: bytes = field(repr=False, default=b"")

    @property
    def start(self):
        return EPOCH + timedelta(seconds=self.start_epoch)

    @property
    def end(self):
        return self.start + timedelta(microseconds=self.n_samples * self.dt_us)


def read_decoded(g1_path):
    g1_path = Path(g1_path)
    stem = g1_path.name[:-len(DECODED_EXTENSIONS["g1"])]
    files = {c: g1_path.with_name(stem + ext) for c, ext in DECODED_EXTENSIONS.items()}
    with open(g1_path, "rb") as f:
        h = f.read(DECODED_HEADER_BYTES)
    if len(h) != DECODED_HEADER_BYTES or _u32(h, 0) != 0xFFFFFFFF or h[6] != 0x80:
        raise ValueError(f"{g1_path.name}: not a decoded EDR file (no start-log header)")
    n = min((p.stat().st_size - DECODED_HEADER_BYTES) // 4 for p in files.values())
    return DecodedRun(folder=g1_path.parent, stem=stem, files=files, start_epoch=_u32(h, 8), stop_epoch=_u32(h, 12),
                      serial=_i32(h, 60),
                      dt_us=struct.unpack(">H", h[68:70])[0], sync0=_u32(h, 206), sync1=_u32(h, 210),
                      skew_us=_i32(h, 214), n_samples=n, header=h)


def _trunc_div(a, b):
    """Integer division rounding toward zero (C behaviour)."""
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b > 0) else -q


def drift(run, t_us, skew_us=None):
    """(shift in samples, rest in us) of the drift correction for a trace starting at t_us (us since 1970)."""
    skew = run.skew_us if skew_us is None else skew_us
    total_s = run.sync1 - run.sync0
    at_start = _trunc_div(skew * (run.start_epoch - run.sync0), total_s)
    rate = float(np.float32(skew) / np.float32(total_s))
    skew_t = math.trunc(at_start + rate * (t_us - run.start_epoch * 1_000_000) / 1_000_000)
    shift = _trunc_div(skew_t, run.dt_us)
    return shift, skew_t - shift * run.dt_us


def text_header(lines):
    """3200-byte EBCDIC textual header from 40 lines of text (e.g. header_edr.txt)."""
    lines = [l.rstrip("\r\n") for l in lines][:40]
    lines += [""] * (40 - len(lines))
    return "".join(l[:80].ljust(80) for l in lines).encode("cp500")


def binary_header(n_samples, dt_us):
    h = bytearray(400)
    for byte_1_based, value in BINARY_HEADER_WORDS.items():
        h[byte_1_based - 3201:byte_1_based - 3199] = struct.pack(">H", value)
    for byte_1_based in (3217, 3219):
        h[byte_1_based - 3201:byte_1_based - 3199] = struct.pack(">H", dt_us)
    for byte_1_based in (3221, 3223):
        h[byte_1_based - 3201:byte_1_based - 3199] = struct.pack(">H", n_samples | 0x8000)
    return bytes(h)


@dataclass
class Receiver:
    """The values EDR takes from Trilobit_RCV (first RC line)."""
    station: int
    line: int = 1
    x: float = 0.0
    y: float = 0.0
    skew_us: int | None = None        # None = skew from the decoded header (= get_skew.ps1)


def trace_header(run, comp, rcv, number, t_us, n_samples, rest_us):
    t = EPOCH + timedelta(microseconds=t_us)
    h = bytearray(240)
    h[0:4] = h[4:8] = h[16:20] = struct.pack(">i", number)
    h[8:12] = struct.pack(">i", rcv.line)
    h[24:28] = struct.pack(">i", run.skew_us if rcv.skew_us is None else rcv.skew_us)
    h[28:30] = struct.pack(">H", CHANNEL_CODE[comp])
    h[52:56] = struct.pack(">I", run.sync0)
    h[56:60] = struct.pack(">I", run.sync1)
    h[70:72] = struct.pack(">h", -100)                       # coordinates in 1/100 m
    h[80:84] = struct.pack(">i", round(rcv.x * 100))
    h[84:88] = struct.pack(">i", round(rcv.y * 100))
    h[88:90] = struct.pack(">H", 1)
    h[90:92] = struct.pack(">H", run.serial)
    h[92:94] = struct.pack(">H", rcv.station)
    h[114:116] = struct.pack(">H", n_samples)
    h[116:118] = struct.pack(">H", run.dt_us)
    # bytes 125-130: accelerometer Gx/Gy/Gz [g * 10000] (header_edr.txt C26-C27); how EDR derives them
    # from aux5 is not known yet -> 0
    h[156:158] = struct.pack(">H", t.year)
    h[158:160] = struct.pack(">H", int(t.strftime("%j")))
    h[160:166] = struct.pack(">HHH", t.hour, t.minute, t.second)
    h[166:168] = struct.pack(">H", 2)                        # UTC
    h[168:172] = struct.pack(">HH", t.microsecond // 1000, t.microsecond % 1000)
    h[176:178] = struct.pack(">h", rest_us)
    return bytes(h)


def trace_samples(run, data, t_us, n_samples, shift):
    """Samples of one trace starting at t_us: decoded data shifted by `shift` samples, zeros outside."""
    first = (t_us - run.start_epoch * 1_000_000) // run.dt_us + shift
    out = np.zeros(n_samples, dtype=">f4")
    a, b = max(first, 0), min(first + n_samples, len(data))
    if b > a:
        out[a - first:b - first] = data[a:b]
    return out


def write_segy(run, comp, out_path, rcv, first_shot_us, n_shots, trace_us, header_lines):
    """One SEG-Y file like EDR's <stem>.<comp>_segy: n_shots traces starting at first_shot_us."""
    n_samples = trace_us // run.dt_us
    data = np.memmap(run.files[comp], dtype=">i4", mode="r", offset=DECODED_HEADER_BYTES, shape=(run.n_samples,))
    with open(out_path, "wb") as f:
        f.write(text_header(header_lines))
        f.write(binary_header(n_samples, run.dt_us))
        for i in range(n_shots):
            t_us = first_shot_us + i * trace_us
            shift, rest = drift(run, t_us, rcv.skew_us)
            f.write(trace_header(run, comp, rcv, i + 1, t_us, n_samples, rest))
            f.write(trace_samples(run, data, t_us, n_samples, shift).tobytes())
    return out_path


# ------------------------------------------------------------------ auxiliary channels
# "Decode Raw" also writes <stem>.aux0/1/2/5: the same 512-byte header, then records of
# epoch [u32, s] + n int16 values. 2026 settings: aux0 every 600 s (looks like temperature in °C:
# 0-2 on the seafloor, 13-19 on deck), aux1 every 60 s (battery voltage [mV]),
# aux2 every 600 s (meaning unknown), aux5 every 1 s (accelerometer x/y/z, 12 bit, ~1024 counts = 1 g).
AUX_VALUES = {"aux0": 1, "aux1": 1, "aux2": 1, "aux5": 3}


def read_aux(path, n_values):
    """(epoch [s] int64 array, int16 values n x n_values) of one aux file."""
    b = np.fromfile(path, dtype=np.uint8, offset=DECODED_HEADER_BYTES)
    size = 4 + 2 * n_values
    b = b[:len(b) // size * size].reshape(-1, size)
    return b[:, :4].copy().view(">u4").ravel().astype(np.int64), b[:, 4:].copy().view(">i2").reshape(-1, n_values)


def aux_per_minute(runs):
    """Aux channels of decoded exports per minute: dict of columns, or None without aux5.

    acc_x/y/z = median of the one-second accelerometer values of the minute (robust against the
    single saturated values at +-2048), tilt_deg = angle between the z axis and the vertical,
    aux0/aux1/aux2 = last value at or before the end of the minute."""
    channels = {}
    for name, n in AUX_VALUES.items():
        parts = [read_aux(p, n) for r in runs if (p := r.folder / f"{r.stem}.{name}").exists()]
        if parts:
            ep = np.concatenate([e for e, _ in parts])
            v = np.concatenate([x for _, x in parts])
            ep, first = np.unique(ep, return_index=True)        # overlapping chunks: each second once
            channels[name] = (ep, v[first])
    if "aux5" not in channels:
        return None
    ep, acc = channels["aux5"]
    minute = ep // 60
    keys, starts = np.unique(minute, return_index=True)
    med = np.array([np.median(g, axis=0) for g in np.split(acc.astype(np.float64), starts[1:])])
    cols = {"time_utc": [EPOCH + timedelta(seconds=int(k) * 60) for k in keys],
            "acc_x": med[:, 0], "acc_y": med[:, 1], "acc_z": med[:, 2],
            "tilt_deg": np.degrees(np.arctan2(np.hypot(med[:, 0], med[:, 1]), med[:, 2]))}
    for name, label in (("aux1", "battery_mV"), ("aux0", "aux0"), ("aux2", "aux2")):
        if name in channels:
            e, v = channels[name]
            i = np.searchsorted(e, (keys + 1) * 60, "left") - 1     # last sample before the minute ends
            cols[label] = np.where(i >= 0, v[np.maximum(i, 0), 0], np.nan)
    return cols
