"""Text report, written into the input folder."""
from datetime import datetime
from pathlib import Path

from . import __version__, core


def write_report(root, scans, settings, results):
    root = Path(root)
    path = root / f"trilobit_report_{datetime.now():%Y%m%d_%H%M%S}.txt"
    lines = [
        f"Trilobit processor {__version__} - {datetime.now():%Y-%m-%d %H:%M:%S}",
        f"Input: {root}",
        f"Segment length: {settings.segment_hours:g} h -> folder '{core.segment_label(settings.segment_hours)}', "
        f"output traces: {settings.output_trace_seconds:g} s",
        f"Start (UTC): {settings.start_utc or 'start of recording'}, end (UTC): {settings.end_utc or 'end of recording'}",
        "",
        "CHECK (EDR export continuous?)",
        core.check_report(scans),
    ]
    if results is not None:
        lines += ["", "OUTPUT"]
        for r in results:
            gap = f", {r['gap_seconds']:.3f} s without data filled with zeros" if r["gap_seconds"] else ""
            lines.append(f"  {r['station']}: {r['files']} files, {r['traces']} traces per component -> {r['folder']}{gap}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
