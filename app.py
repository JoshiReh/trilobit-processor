"""Trilobit Processor - Mac app (tkinter).

Pick the data folder (drag and drop or "Browse..."), choose the segment length,
press "Check and process". Check, header fix, Galperin rotation and segmentation
run in a background thread; the window only shows progress and the log.
"""
import queue
import subprocess
import threading
import time
import tkinter as tk
from datetime import timedelta
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from trilobit_processor import __version__, core
from trilobit_processor.report import write_report

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAS_DND = True
except Exception:          # app still works without drag and drop
    HAS_DND = False

SEGMENT_CHOICES = ["1", "6", "12", "24", "48", "72"]


class App:
    def __init__(self, root):
        self.root = root
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.worker = None
        self.scans = None

        root.title("Trilobit Processor")
        root.minsize(720, 620)
        main = ttk.Frame(root, padding=14)
        main.pack(fill="both", expand=True)

        # --- input folder
        box = ttk.LabelFrame(main, text="Input folder (EDR exports *.g1/g2/g3/h1_segy anywhere inside, any names)",
                             padding=10)
        box.pack(fill="x")
        self.drop = tk.Label(box, text="Drop the folder here" if HAS_DND else "Choose the folder with “Browse…”",
                             relief="groove", bd=2, height=3, fg="gray40")
        self.drop.pack(fill="x", pady=(0, 8))
        row = ttk.Frame(box)
        row.pack(fill="x")
        self.path = tk.StringVar()
        ttk.Entry(row, textvariable=self.path).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Browse…", command=self.browse).pack(side="left", padx=(8, 0))
        if HAS_DND:
            self.drop.drop_target_register(DND_FILES)
            self.drop.dnd_bind("<<Drop>>", self.on_drop)

        # --- settings
        box = ttk.LabelFrame(main, text="Settings", padding=10)
        box.pack(fill="x", pady=10)
        self.hours = tk.StringVar(value="24")
        self.trace_s = tk.StringVar(value="60")
        self.start = tk.StringVar()
        self.end = tk.StringVar()
        self.folder_hint = tk.StringVar()

        ttk.Label(box, text="Segment length [hours]").grid(row=0, column=0, sticky="w")
        ttk.Combobox(box, textvariable=self.hours, values=SEGMENT_CHOICES, width=8).grid(row=0, column=1, sticky="w")
        ttk.Label(box, textvariable=self.folder_hint, foreground="gray40").grid(row=0, column=2, columnspan=3, sticky="w", padx=10)
        ttk.Label(box, text="Trace length [s]").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(box, textvariable=self.trace_s, width=8).grid(row=1, column=1, sticky="w")
        ttk.Label(box, text="Start (optional, UTC)").grid(row=2, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.start, width=20).grid(row=2, column=1, sticky="w")
        ttk.Label(box, text="End (optional, UTC)").grid(row=2, column=2, sticky="w", padx=(16, 0))
        ttk.Entry(box, textvariable=self.end, width=20).grid(row=2, column=3, sticky="w")
        ttk.Label(box, text="Format: YYYY-MM-DD HH:MM[:SS] in UTC, e.g. 2026-07-07 16:00 - empty = whole recording "
                            "(with decoded files: incl. transport and deployment)", foreground="gray40").grid(
            row=3, column=0, columnspan=4, sticky="w")
        self.hours.trace_add("write", lambda *_: self.update_hint())
        self.update_hint()

        # --- buttons + progress
        row = ttk.Frame(main)
        row.pack(fill="x")
        self.btn_check = ttk.Button(row, text="Check only", command=lambda: self.start_check(process=False))
        self.btn_run = ttk.Button(row, text="Check and process", command=lambda: self.start_check(process=True))
        self.btn_cancel = ttk.Button(row, text="Cancel", command=self.cancel.set, state="disabled")
        for b in (self.btn_check, self.btn_run, self.btn_cancel):
            b.pack(side="left", padx=(0, 8))
        self.progress = ttk.Progressbar(main, mode="determinate", maximum=1.0)
        self.progress.pack(fill="x", pady=(10, 2))
        self.status = tk.StringVar(value="Ready.")
        ttk.Label(main, textvariable=self.status).pack(anchor="w")

        self.log_box = ScrolledText(main, height=16, font=("Menlo", 11), state="disabled")
        self.log_box.pack(fill="both", expand=True, pady=(8, 0))
        self.log(f"Trilobit Processor {__version__}")
        root.after(100, self.poll)

    # ---------------------------------------------------------------- helpers
    def log(self, text):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", text + "\n")
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    def update_hint(self):
        try:
            self.folder_hint.set(f"→ output folder per station: “{core.segment_label(float(self.hours.get()))}”")
        except ValueError:
            self.folder_hint.set("")

    def browse(self):
        p = filedialog.askdirectory(title="Choose the folder with the station folders")
        if p:
            self.path.set(p)

    def on_drop(self, event):
        paths = self.root.tk.splitlist(event.data)
        if paths:
            p = Path(paths[0])
            self.path.set(str(p if p.is_dir() else p.parent))

    def busy(self, running):
        state = "disabled" if running else "normal"
        self.btn_check.configure(state=state)
        self.btn_run.configure(state=state)
        self.btn_cancel.configure(state="normal" if running else "disabled")

    def read_settings(self):
        hours = float(self.hours.get())
        trace_s = float(self.trace_s.get())
        if hours <= 0 or trace_s <= 0:
            raise ValueError("Segment and trace length must be greater than 0.")
        start, end = core.parse_utc(self.start.get()), core.parse_utc(self.end.get())
        if start and end and end <= start:
            raise ValueError("The end must be after the start.")
        return core.Settings(hours, trace_s, start, end)

    def run_thread(self, target, *args):
        self.cancel.clear()
        self.busy(True)
        self.started = time.monotonic()
        self.worker = threading.Thread(target=target, args=args, daemon=True)
        self.worker.start()

    def progress_cb(self, frac, text):
        self.events.put(("progress", frac, text))

    def log_cb(self, text):
        self.events.put(("log", text))

    def eta(self, frac):
        """' - 0:04:12 elapsed, about 0:20:00 left' (left only once there is enough progress)."""
        elapsed = time.monotonic() - self.started
        fmt = lambda s: str(timedelta(seconds=int(s)))
        text = f" - {fmt(elapsed)} elapsed"
        if frac > 0.01:
            text += f", about {fmt(elapsed * (1 - frac) / frac)} left"
        return text

    # ---------------------------------------------------------------- step 0
    def start_check(self, process):
        root = Path(self.path.get().strip())
        if not root.is_dir():
            messagebox.showerror("No folder", "Please choose a valid input folder first.")
            return
        try:
            settings = self.read_settings()
        except ValueError as e:
            messagebox.showerror("Settings", f"Invalid input: {e}")
            return
        stations = core.find_stations(root)
        if not stations:
            messagebox.showerror("No EDR exports", "This folder contains no EDR exports (*.g1_segy files).")
            return
        self.log(f"\n=== Checking {len(stations)} station(s) in {root} ===")
        self.run_thread(self.check_worker, root, settings, process)

    def check_worker(self, root, settings, process):
        try:
            scans = core.scan_all(root, progress=self.progress_cb, log=self.log_cb)
            self.events.put(("checked", root, settings, scans, process))
        except Exception as e:
            self.events.put(("error", f"Error during the check: {e}"))

    def after_check(self, root, settings, scans, process):
        self.log(core.check_report(scans))
        ok = all(s.ok for s in scans)
        if not process:
            path = write_report(root, scans, settings, None)
            self.finish("Check OK." if ok else "Check: PROBLEMS found.", path)
            return
        if not ok:
            go = messagebox.askyesno(
                "Problems found",
                "The check found gaps or mismatches (see log).\n\n"
                "Recommended: export again in EDR and check again.\n\n"
                "Process anyway? Missing times are filled with zeros, "
                "all data stay at their correct time.", icon="warning")
            if not go:
                path = write_report(root, scans, settings, None)
                self.finish("Stopped after the check.", path)
                return
        existing = core.existing_outputs(scans, settings)
        if existing:
            folders = sorted({str(f.parent) for f in existing})
            if not messagebox.askyesno("Output exists",
                                       f"{len(existing)} segment files already exist in:\n\n"
                                       + "\n".join(folders[:8]) + ("\n…" if len(folders) > 8 else "")
                                       + "\n\nOverwrite?"):
                self.finish("Stopped, nothing overwritten.", None)
                return
        plans = [p for p in (core.plan_station(s, settings) for s in scans) if p is not None]
        self.log("\n=== Overview ===")
        self.log(core.processing_summary(plans, settings))
        self.cancel.clear()
        self.busy(True)
        self.status.set("Starting in 10 s …")
        self.root.after(10_000, lambda: self.start_processing(root, settings, scans, existing, plans))

    def start_processing(self, root, settings, scans, existing, plans):
        if self.cancel.is_set():
            self.finish("Cancelled, nothing written.", None)
            return
        self.log("\n" + core.processing_sentence(plans, settings))
        self.log("\n=== Processing ===")
        self.run_thread(self.process_worker, root, settings, plans, existing)

    # ---------------------------------------------------------------- steps 1+2
    def process_worker(self, root, settings, plans, existing):
        try:
            for f in existing:
                f.unlink()
            scans = [p.scan for p in plans]
            total = sum(p.n_output_traces() for p in plans) or 1
            done = 0
            results = []
            for plan in plans:
                n = plan.n_output_traces()
                cb = lambda frac, text, d=done, n=n: self.progress_cb((d + frac * n) / total, text)
                results.append(core.process_station(plan, settings, progress=cb, cancel=self.cancel, log=self.log_cb))
                done += n
                r = results[-1]
                gap = f", {r['gap_seconds']:.3f} s without data (zeros)" if r["gap_seconds"] else ""
                self.events.put(("log", f"   → {r['files']} files in {r['folder']}{gap}"))
            path = write_report(root, scans, settings, results)
            self.events.put(("done", "Done.", path))
        except core.Cancelled:
            self.events.put(("done", "Cancelled (the unfinished segment was deleted).", None))
        except Exception as e:
            self.events.put(("error", f"Error during processing: {e}"))

    # ---------------------------------------------------------------- event loop
    def finish(self, text, report_path):
        self.busy(False)
        self.status.set(text)
        self.log(text)
        if report_path:
            self.log(f"Report: {report_path}")
            self.progress["value"] = 1.0

    def poll(self):
        try:
            while True:
                ev = self.events.get_nowait()
                kind = ev[0]
                if kind == "progress":
                    self.progress["value"] = ev[1]
                    self.status.set(ev[2] + self.eta(ev[1]))
                elif kind == "log":
                    self.log(ev[1])
                elif kind == "checked":
                    self.busy(False)
                    self.after_check(*ev[1:])
                elif kind == "done":
                    self.finish(ev[1], ev[2])
                    if ev[2]:
                        subprocess.run(["open", "-R", str(ev[2])], check=False)  # show report in Finder
                elif kind == "error":
                    self.finish(ev[1], None)
                    messagebox.showerror("Error", ev[1])
        except queue.Empty:
            pass
        self.root.after(100, self.poll)


def main():
    global HAS_DND
    try:
        root = TkinterDnD.Tk() if HAS_DND else tk.Tk()
    except Exception:      # tkdnd library could not be loaded
        HAS_DND = False
        root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
