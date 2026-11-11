# <img src="assets/icon.png" width="48" align="center" alt=""> Trilobit Processor

Turns the decoded data of a Trilobit ocean bottom node (EDRTest208 "Decode Raw") into SEG-Y files of any
length - in one step, without the rest of EDRTest208.

<!-- screen recording (20 s) goes here, e.g. ![How it works](docs/demo.gif) -->

## Download

| | |
|---|---|
| **macOS** (Apple Silicon, M1 or newer) | `TrilobitProcessor-mac.zip` under [Releases](../../releases/latest) |
| **Windows** (10/11) | `TrilobitProcessor-windows.zip` under [Releases](../../releases/latest) |

Unzip and start the app. The first time, your computer may block it because it is not signed:

- **macOS:** right-click the app → *Open* → *Open*.
- **Windows:** unzip the whole folder, start `Trilobit Processor.exe`, and if a blue window
  appears click *More info* → *Run anyway*.

## How to use

1. In EDRTest208: download the data, run `fix_all` and **Decode Raw** (whole recording, no chunking).
2. Drop the folder with the decoded files (`.g1 .g2 .g3 .h1`, `.aux*`, `.txt`) into the app.
   The folder name must contain the 4-digit station number, e.g. `full0989`.
3. Choose the segment length (e.g. 24 h) and, optionally, a start and end time in UTC.
4. Click **Check and process**.

The app checks the data, corrects the clock drift exactly like EDRTest208, and writes:

- `24hourChunks/` - one SEG-Y file per segment and component (x, y, z, h)
- `<station>_aux_per_minute.csv` - tilt, battery and temperature for every minute
- `trilobit_report_<date>.txt` - what was checked and written

## Jupyter notebooks

If you exported SEG-Y files with EDRTest208 yourself, the notebooks in [notebooks/](notebooks/) do the same
steps one by one (export check, header fix and x/y/z, daily files, quality check). They are used in Session 6 of
the [Geo-OBS manual](https://github.com/JoshiReh/Geo-OBS).

## More

Technical details, command line version, tests and how to build: [docs/TECHNICAL.md](docs/TECHNICAL.md)

## Credits

Written by Joshua Rehlinghaus. The file formats and the drift correction are based on the EDRTest208 tools
and helper scripts by **Bent Ole Ruud** (University of Bergen), e.g. `fix_rollover` and `get_skew`.

## License

MIT, see [LICENSE](LICENSE).
