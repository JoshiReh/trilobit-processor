"""Command line version: python -m trilobit_processor.cli /path/to/data --hours 24"""
import argparse
import sys

from . import core
from .report import write_report


def main(argv=None):
    p = argparse.ArgumentParser(description="Check, fix and segment Trilobit EDR exports.")
    p.add_argument("root", help="folder with EDR exports (*.g1_segy ...) or decoded files (*.g1 ...) anywhere inside")
    p.add_argument("--hours", type=float, default=24, help="segment length in hours (default 24)")
    p.add_argument("--trace-seconds", type=float, default=60, help="output trace length in s (default 60)")
    p.add_argument("--start", help="optional start, UTC 'YYYY-MM-DD HH:MM[:SS]'")
    p.add_argument("--end", help="optional end, UTC 'YYYY-MM-DD HH:MM[:SS]'")
    p.add_argument("--quiet", action="store_true", help="no step-by-step log")
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--force", action="store_true", help="process even if the check finds problems")
    p.add_argument("--overwrite", action="store_true", help="replace existing segment files")
    a = p.parse_args(argv)

    settings = core.Settings(a.hours, a.trace_seconds, core.parse_utc(a.start), core.parse_utc(a.end))
    log = (lambda text: None) if a.quiet else print

    scans = core.scan_all(a.root, log=log)
    print(core.check_report(scans))
    if a.check_only:
        print("report:", write_report(a.root, scans, settings, None))
        return 0 if all(s.ok for s in scans) else 1
    if not all(s.ok for s in scans) and not a.force:
        print("\nStopped: the check found problems (use --force to process anyway, gaps are filled with zeros).")
        print("report:", write_report(a.root, scans, settings, None))
        return 1

    existing = core.existing_outputs(scans, settings)
    if existing and not a.overwrite:
        print(f"\nStopped: {len(existing)} segment files already exist (use --overwrite).")
        return 1
    for f in existing:
        f.unlink()

    plans = [p for p in (core.plan_station(s, settings) for s in scans) if p is not None]
    print(core.processing_summary(plans, settings))
    print(core.processing_sentence(plans, settings))
    results = [core.process_station(plan, settings, log=log) for plan in plans]
    print("report:", write_report(a.root, scans, settings, results))
    return 0


if __name__ == "__main__":
    sys.exit(main())
