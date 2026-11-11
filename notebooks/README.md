# Notebooks: the classic route (Geo-OBS manual, Session 6)

Python notebooks for processing Trilobit ocean bottom seismometer data, used in Session 6 of the Geo-OBS manual
(https://github.com/JoshiReh/Geo-OBS). They do the same as the app, step by step and for SEG-Y files
exported with EDRTest208 ("Seis to SEGY"). The app and the notebooks give identical data. Run them in this order:

| Step | Notebook | What it does |
|---|---|---|
| 1 | `1_check_edr_export.ipynb` | Checks the EDR export for gaps and wrong clock skew |
| 2 | `2_fix_header_and_xyz.ipynb` | Repairs the binary header and converts g1/g2/g3 to x/y/z |
| 3 | `3_cut_daily_files.ipynb` | Joins the chunks and cuts them into daily files |
| 4 | `4_qc_daily_files.ipynb` | Plots one day to check the data quality |

Only the cell **Settings** at the top of each notebook has to be changed.
Requirements: Python 3 with numpy, scipy, pandas, matplotlib and jupyter.
