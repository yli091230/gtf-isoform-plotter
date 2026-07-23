"""Allow ``python -m transcript_isoform_plotter`` execution."""

from .cli import main


if __name__ == "__main__":
    raise SystemExit(main())

