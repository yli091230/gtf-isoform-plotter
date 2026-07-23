#!/usr/bin/env python3
"""Backward-compatible launcher for gtf-isoform-plotter."""

from transcript_isoform_plotter.cli import main


if __name__ == "__main__":
    raise SystemExit(main())

