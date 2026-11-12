"""Verify a delivered offline report without reading raw expression or running inference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .exporter import _sha256, validate_bundle


def verify_report(directory: Path) -> dict[str, Any]:
    directory = directory.resolve()
    manifest = json.loads((directory / "REPORT_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("statistics_recomputed") is not False:
        raise ValueError("Report manifest does not attest to renderer purity")
    for item in manifest["outputs"]:
        path = (directory / item["path"]).resolve()
        if not path.is_relative_to(directory):
            raise ValueError("Report manifest path escapes delivery directory")
        if path.stat().st_size != item["bytes"] or _sha256(path) != item["sha256"]:
            raise ValueError(f"Report output checksum mismatch: {item['path']}")
    bundle = json.loads((directory / "downloads/bundle.json").read_text(encoding="utf-8"))
    validate_bundle(bundle)
    counts = {key: len(bundle[key]) for key in
              ["units", "genetic_links", "donor_evidence", "sections"]}
    if counts != manifest["counts"]:
        raise ValueError("Report count mismatch")
    return {"status": "PASS_REPORT_DELIVERY_INTEGRITY_ONLY", "counts": counts,
            "statistics_recomputed": False, "files_verified": len(manifest["outputs"]),
            "browser_visual_acceptance": "NOT_ESTABLISHED_BY_THIS_CHECK",
            "scientific_calibration": "NOT_ESTABLISHED_BY_THIS_CHECK"}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify_report(args.directory), indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
