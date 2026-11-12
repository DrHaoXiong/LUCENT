"""Export an existing report bundle: python -m lucent.report --help."""

from .exporter import main

if __name__ == "__main__":
    raise SystemExit(main())
