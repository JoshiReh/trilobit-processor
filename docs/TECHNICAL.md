# Trilobit Processor - technical details

Desktop app (macOS, Windows) that turns Trilobit recordings into time segments, in one go. It takes either the output of
EDR "Decode Raw" (recommended: then EDR is only needed for the download, `fix_all_ln.ps1` and
"Decode Raw") or SEG-Y files from EDR "Seis to SEGY":

0. **Convert** (decoded files only) - does what EDR "Seis to SEGY" does: 60 s traces, clock drift
   correction from the skew in the file header, EDR's headers. Checked against EDR byte for byte
   (see `trilobit_processor/edr.py`); no `Trilobit_SHT`/`Trilobit_RCV`/`get_skew` needed
1. **Check** - is every station's EDR export one continuous time series? (no gaps inside a file,
   no gaps between exports, start/stop = EDR log, g1/g2/g3/h1 identical in time)
2. **Fix** - repairs the broken EDR binary header and rotates the Galperin components g1/g2/g3 into x/y/z
3. **Segment** - cuts the continuous recording into segments of any length (1 h, 12 h, 24 h, 72 h, ...)
   with 60 s traces

Every trace is placed by the time stamp in its own header (a gap stays a gap, filled with zeros,
instead of shifting all later data), and the input folders are only read - nothing is changed or
deleted there.

## Input

Any folder that contains exports, at any depth, with any folder and file names. One export = four
files with the same name and the endings

- `.g1`, `.g2`, `.g3`, `.h1` - output of EDR "Decode Raw" (unchunked, the whole recording), or
- `.g1_segy`, `.g2_segy`, `.g3_segy`, `.h1_segy` - output of EDR "Seis to SEGY".

Optional: `<same name>.txt` = EDR log, and a `.raw` file in the same folder or the station folder,
both only used by the check. A station can have one export or many (chunks). If both kinds exist
with the same name, the decoded files are used.

Decoded files contain the whole recording, including transport and deployment: set **start and end
(UTC)** in the app. The first trace starts exactly at the start (like `Trilobit_SHT` in EDR).

```
<input folder>/                       e.g. several stations ...
    0498/
        1stChunk/
            0498_1stChunk_seq1.g1_segy   (.g2_segy, .g3_segy, .h1_segy, .txt)
        2ndChunk/ ...
    0578/ ...

adbasdbad_0989/                       ... or a single, unchunked station
    full0989.g1                       (.g2, .g3, .h1, .aux*, .txt from "Decode Raw")
```

Decoded files have no station number: it is taken from the names, the 4-digit number closest to the
file (e.g. `full0989.g1` or `0989/whatever.g1`). For SEG-Y exports the station number is read from
the trace headers (EDR writes it from `Trilobit_RCV`), so exports are grouped by station whatever
their names. If a folder or file name contains a 4-digit number
(e.g. `0989`), the check reports a problem when it does not match the station number in the
headers (a `Trilobit_RCV` copied from another station).

## Output

Per station a new folder in the deepest folder that contains all of its exports (e.g. `0498/`):
`<hours>hourChunks`, e.g. `24hourChunks`, `12hourChunks`, `1hourChunks`.
Files: `<station>_<x|y|z|h>_chunk<NNNN>_<start UTC>.segy`.
A report `trilobit_report_<date>.txt` is written into the input folder.

With decoded files there is also `<station>_aux_per_minute.csv` in the station folder, for the whole
recording (not cut to start/end), from the aux files of "Decode Raw":

| column | meaning |
|---|---|
| `time_utc` | start of the minute |
| `acc_x`, `acc_y`, `acc_z` | accelerometer (`aux5`, 12 bit, ~1024 counts = 1 g), median of the 60 one-second values |
| `tilt_deg` | angle between the node's z axis and the vertical, from acc_x/y/z (independent of the scale) |
| `battery_mV` | battery voltage (`aux1`) |
| `aux0`, `aux2` | raw values. aux0 looks like the temperature in °C (0989: 13-15 on deck, 0-2 on the seafloor, 17-19 after recovery); aux2 unknown |

EDR writes accelerometer values into its SEG-Y trace headers too (bytes 125-130, "ACC Gx/Gy/Gz (G*10000)"
in header_edr.txt); how it derives them from aux5 is not known, the app writes 0 there.

If the check finds a problem, the app asks before processing. Recommended: export again in EDR.

## Use

- App: drag the input folder into the window (or "Browse..."), choose the segment length,
  "Check and process".
- Command line: `python -m trilobit_processor.cli /Volumes/.../zJoshi --hours 12`
  (`--check-only`, `--force`, `--overwrite`, `--start/--end "YYYY-MM-DD HH:MM[:SS]"` in UTC, `--quiet`)
- The log shows every step (files found, check per export, every segment written) and the status line
  the elapsed and remaining time.

## Build

Releases are built by GitHub Actions (`.github/workflows/build.yml`, macOS and Windows) when a
version tag is pushed. Locally on a Mac:

```
conda create -n trilobit-app python=3.11 numpy tk pyinstaller
conda run -n trilobit-app pip install tkinterdnd2
./build_app.sh            # -> dist/Trilobit Processor.app (macOS)
```

## Tests

```
python tests/test_pipeline.py /tmp/trilobit_tests [/Volumes/NilRodes/zJoshi/1069]
```

Synthetic data with the 2026 problem (59 s missing between two exports), wrong `Trilobit_RCV` values,
any names, decoded files (decoded -> app must equal EDR SEG-Y -> app), and optionally a real
station: its first 3 h are turned back into raw EDR-like files and the output is compared with the
existing output of the same station (`DailyChunk`).

```
python tests/verify_edr.py <chunk folder> ...
```

Compares `edr.py` with EDR's own SEG-Y, trace by trace, in folders that contain the decoded files,
EDR's `*_segy`, `Trilobit_SHT`, `Trilobit_RCV` and `header_edr.txt` (2026: identical for 0676
1st/2nd/3rdChunk and 0989 2ndChunk, except header bytes 125-130 and EDR reading its own file header
as samples when the drift correction points before the first sample).
