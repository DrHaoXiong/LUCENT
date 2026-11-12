"""Validate and serialize existing LUCENT evidence into a portable offline report.

Only presentation data are assembled here. No p values, adjustments, spatial
edges, expression values, or scores are estimated by this module. The package
deliberately uses the Python standard library only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class ReportValidationError(ValueError):
    """The input cannot be represented faithfully by the report contract."""


_ASSETS = Path(__file__).parent / "assets"
_SHA256 = re.compile(r"[0-9a-fA-F]{64}\Z")


def _fail(where: str, message: str) -> None:
    raise ReportValidationError(f"{where}: {message}")


def _object(value: Any, where: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        _fail(where, "expected an object")
    return value


def _array(value: Any, where: str) -> list[Any]:
    if not isinstance(value, list):
        _fail(where, "expected an array")
    return value


def _text(value: Any, where: str, *, nullable: bool = False) -> str | None:
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not value.strip():
        _fail(where, "expected a nonempty string")
    return value


def _number(value: Any, where: str, *, nullable: bool = True,
            probability: bool = False) -> int | float | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(where, "expected a finite number" + (" or null" if nullable else ""))
    try:
        finite = math.isfinite(value)
    except OverflowError:
        finite = False
    if not finite:
        _fail(where, "number must be finite")
    if probability and not 0 <= value <= 1:
        _fail(where, "probability must lie in [0, 1]")
    return value


def _count(value: Any, where: str, *, nullable: bool = False) -> None:
    if value is None and nullable:
        return
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        _fail(where, "expected a nonnegative integer")


def _required(row: dict[str, Any], names: tuple[str, ...], where: str) -> None:
    for name in names:
        if name not in row:
            _fail(where, f"missing {name!r}")


def _unique(value: Any, seen: set[Any], where: str) -> None:
    if value in seen:
        _fail(where, f"duplicate identifier or record {value!r}")
    seen.add(value)


def _json_tree(value: Any, where: str = "bundle") -> None:
    """Reject values JSON would silently coerce or JavaScript cannot represent."""
    if value is None or isinstance(value, (str, bool)):
        return
    if isinstance(value, (int, float)):
        _number(value, where, nullable=False)
        if isinstance(value, int) and abs(value) > 9007199254740991:
            _fail(where, "integer exceeds JavaScript's exact integer range")
        return
    if isinstance(value, list):
        for i, item in enumerate(value):
            _json_tree(item, f"{where}[{i}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                _fail(where, "object keys must be strings")
            _json_tree(item, f"{where}.{key}")
        return
    _fail(where, f"not a JSON value: {type(value).__name__}")


def _result(row: Any, where: str, *, has_q: bool) -> None:
    row = _object(row, where)
    fields = ("score", "p", "status") + (("q", "evaluable") if has_q else ())
    _required(row, fields, where)
    _number(row["score"], where + ".score")
    _number(row["p"], where + ".p", probability=True)
    _text(row["status"], where + ".status")
    if has_q:
        _number(row["q"], where + ".q", probability=True)
        if row["q"] is not None and row["p"] is None:
            # Some frozen families deliberately retain unavailable hypotheses
            # with a correction-only placeholder p=1. Preserve their saved q,
            # but require explicit provenance and keep them unevaluable.
            placeholder = row.get("correction_p")
            if row.get("evaluable") is not False or placeholder is None:
                _fail(where + ".q", "an adjustment without a source p value requires an explicit correction_p placeholder and evaluable=false")
            _number(placeholder, where + ".correction_p", nullable=False, probability=True)
        if row["evaluable"] is not None and not isinstance(row["evaluable"], bool):
            _fail(where + ".evaluable", "expected a boolean or null")


def validate_bundle(bundle: dict[str, Any]) -> None:
    """Check types, probability domains, identifiers, and cross references.

    Free text status values and additional provenance fields are retained as
    supplied. Null and zero are never substituted for one another. Spot IDs
    are local to a section; edges must use that section's integer indices.
    """
    _object(bundle, "bundle")
    _json_tree(bundle)
    _required(bundle, ("metadata", "units", "genetic_links", "donor_evidence", "sections"),
              "bundle")
    meta = _object(bundle["metadata"], "metadata")
    _required(meta, ("schema_version", "title", "trait", "dataset", "claim", "source_status",
                     "default_version", "default_section_id", "versions", "sources", "notes"),
              "metadata")
    if meta["schema_version"] != "1.0":
        _fail("metadata.schema_version", "supported version is '1.0'")
    for key in ("title", "trait", "dataset", "claim", "source_status"):
        _text(meta[key], f"metadata.{key}")
    _text(meta["default_version"], "metadata.default_version")
    _text(meta["default_section_id"], "metadata.default_section_id", nullable=True)
    for i, note in enumerate(_array(meta["notes"], "metadata.notes")):
        _text(note, f"metadata.notes[{i}]")
    versions: set[str] = set()
    for i, version in enumerate(_array(meta["versions"], "metadata.versions")):
        where = f"metadata.versions[{i}]"
        version = _object(version, where)
        _required(version, ("id", "label", "status", "p_label", "q_label",
                            "adjustment_method", "family_size"), where)
        for key in ("id", "label", "status", "p_label", "q_label"):
            _text(version[key], where + "." + key)
        _unique(version["id"], versions, where + ".id")
        _text(version["adjustment_method"], where + ".adjustment_method", nullable=True)
        _count(version["family_size"], where + ".family_size", nullable=True)
        if str(version["adjustment_method"]).lower() == "holm":
            if "fdr" in version["q_label"].lower():
                _fail(where + ".q_label", "Holm is adjusted p, not FDR")
    if not versions or meta["default_version"] not in versions:
        _fail("metadata.default_version", "must reference a declared version")
    for i, source in enumerate(_array(meta["sources"], "metadata.sources")):
        where = f"metadata.sources[{i}]"
        source = _object(source, where)
        _required(source, ("label", "path", "sha256"), where)
        _text(source["label"], where + ".label")
        _text(source["path"], where + ".path")
        if source["sha256"] is not None and (
            not isinstance(source["sha256"], str) or not _SHA256.fullmatch(source["sha256"])
        ):
            _fail(where + ".sha256", "expected a SHA-256 hex digest or null")

    unit_by_id: dict[str, dict[str, Any]] = {}
    unit_ids: set[str] = set()
    channels: set[str] = set()
    hypotheses: set[tuple[str, str, str, str]] = set()
    for i, unit in enumerate(_array(bundle["units"], "units")):
        where = f"units[{i}]"
        unit = _object(unit, where)
        fields = ("unit_id", "channel_id", "lr_pair", "context", "sender_state", "receiver_state")
        _required(unit, fields + ("ligand_genes", "receptor_genes", "versions"), where)
        for key in fields:
            _text(unit[key], where + "." + key)
        _unique(unit["unit_id"], unit_ids, where + ".unit_id")
        unit_by_id[unit["unit_id"]] = unit
        channels.add(unit["channel_id"])
        _unique(tuple(unit[k] for k in ("channel_id", "context", "sender_state", "receiver_state")),
                hypotheses, where)
        for key in ("ligand_genes", "receptor_genes"):
            gene_ids: set[str] = set()
            for j, gene in enumerate(_array(unit[key], where + "." + key)):
                _text(gene, f"{where}.{key}[{j}]")
                _unique(gene, gene_ids, where + "." + key)
        results = _object(unit["versions"], where + ".versions")
        for version_id, result in results.items():
            if version_id not in versions:
                _fail(where + ".versions", f"unknown version {version_id!r}")
            _result(result, where + ".versions." + version_id, has_q=True)
    if meta.get("default_unit_id") is not None and meta["default_unit_id"] not in unit_ids:
        _fail("metadata.default_unit_id", "must reference a declared unit")

    section_by_id: dict[str, dict[str, Any]] = {}
    section_ids: set[str] = set()
    for i, section in enumerate(_array(bundle["sections"], "sections")):
        where = f"sections[{i}]"
        section = _object(section, where)
        _required(section, ("id", "donor_id", "context", "label", "coordinate_system",
                            "spots", "edges", "expression", "expression_label"), where)
        for key in ("id", "donor_id", "context", "label", "expression_label"):
            _text(section[key], where + "." + key)
        _unique(section["id"], section_ids, where + ".id")
        section_by_id[section["id"]] = section
        coord = _object(section["coordinate_system"], where + ".coordinate_system")
        _text(coord.get("unit"), where + ".coordinate_system.unit")
        if coord.get("y_direction") not in ("up", "down"):
            _fail(where + ".coordinate_system.y_direction", "must be 'up' or 'down'")
        spots = _array(section["spots"], where + ".spots")
        spot_ids: set[str] = set()
        for j, spot in enumerate(spots):
            sw = f"{where}.spots[{j}]"
            spot = _object(spot, sw)
            for key in ("id", "state"):
                _text(spot.get(key), sw + "." + key)
            _unique(spot["id"], spot_ids, sw + ".id")
            _number(spot.get("x"), sw + ".x", nullable=False)
            _number(spot.get("y"), sw + ".y", nullable=False)
            if "section_id" in spot and spot["section_id"] != section["id"]:
                _fail(sw + ".section_id", "cross-section spot reference")
        edge_ids: set[tuple[int, int]] = set()
        for j, edge in enumerate(_array(section["edges"], where + ".edges")):
            ew = f"{where}.edges[{j}]"
            if not isinstance(edge, list) or len(edge) != 2:
                _fail(ew, "expected two local integer spot indices")
            for index in edge:
                if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(spots):
                    _fail(ew, "edge index lies outside this section's spots")
            _unique(tuple(edge), edge_ids, ew)
        expression = _object(section["expression"], where + ".expression")
        for gene, values in expression.items():
            _text(gene, where + ".expression gene")
            values = _array(values, where + ".expression." + gene)
            if len(values) != len(spots):
                _fail(where + ".expression." + gene, "one value is required per local spot")
            for j, value in enumerate(values):
                _number(value, f"{where}.expression.{gene}[{j}]")
        if section.get("image") is not None:
            image = _object(section["image"], where + ".image")
            _text(image.get("path"), where + ".image.path")
            bounds = _array(image.get("bounds"), where + ".image.bounds")
            if len(bounds) != 4:
                _fail(where + ".image.bounds", "expected [xmin, ymin, xmax, ymax]")
            for j, value in enumerate(bounds):
                _number(value, f"{where}.image.bounds[{j}]", nullable=False)
            xmin, ymin, xmax, ymax = bounds
            if xmin >= xmax or ymin >= ymax:
                _fail(where + ".image.bounds", "bounds must have positive width and height")
            if any(not (xmin <= s["x"] <= xmax and ymin <= s["y"] <= ymax) for s in spots):
                _fail(where + ".image.bounds", "registered bounds must contain all section coordinates")
    default_section = meta["default_section_id"]
    if default_section is not None and default_section not in section_ids:
        _fail("metadata.default_section_id", "must reference a declared section or be null")
    if section_ids and default_section is None:
        _fail("metadata.default_section_id", "a default section is required when sections exist")

    genetic_ids: set[tuple[Any, ...]] = set()
    for i, link in enumerate(_array(bundle["genetic_links"], "genetic_links")):
        where = f"genetic_links[{i}]"
        link = _object(link, where)
        _text(link.get("channel_id"), where + ".channel_id")
        if link["channel_id"] not in channels:
            _fail(where + ".channel_id", "unknown channel")
        for key in ("locus_id", "variant_id", "rsid", "gene_id", "gene_symbol", "endpoint_role",
                    "mapping_basis", "pip_status", "credible_set", "build"):
            if key in link:
                _text(link[key], where + "." + key, nullable=True)
        for key in ("pip", "gwas_p"):
            _number(link.get(key), where + "." + key, probability=True)
        _unique(tuple(link.get(k) for k in ("channel_id", "locus_id", "variant_id", "rsid",
                                           "gene_id", "gene_symbol", "endpoint_role", "mapping_basis",
                                           "credible_set", "build")), genetic_ids, where)
        if link.get("pip") is None and not link.get("pip_status"):
            _fail(where + ".pip_status", "a null PIP requires an explicit status")

    donor_ids: set[tuple[str, str, str]] = set()
    for i, evidence in enumerate(_array(bundle["donor_evidence"], "donor_evidence")):
        where = f"donor_evidence[{i}]"
        evidence = _object(evidence, where)
        _required(evidence, ("unit_id", "version", "donor_id", "sections", "edges", "loo"), where)
        for key in ("unit_id", "version", "donor_id"):
            _text(evidence[key], where + "." + key)
        if evidence["unit_id"] not in unit_ids:
            _fail(where + ".unit_id", "unknown unit")
        if evidence["version"] not in unit_by_id[evidence["unit_id"]]["versions"]:
            _fail(where + ".version", "version does not exist on referenced unit")
        _unique(tuple(evidence[k] for k in ("unit_id", "version", "donor_id")), donor_ids, where)
        _result(evidence, where, has_q=False)
        unit_context = unit_by_id[evidence["unit_id"]]["context"]
        matching_sections = {s["id"] for s in section_by_id.values()
                             if s["donor_id"] == evidence["donor_id"] and s["context"] == unit_context}
        if evidence["donor_id"] not in {s["donor_id"] for s in section_by_id.values()}:
            _fail(where + ".donor_id", "unknown donor")
        if not matching_sections and (evidence["sections"] not in (0, []) or evidence["edges"] != 0):
            _fail(where + ".donor_id", "nonzero evidence has no section for this donor and unit context")
        if isinstance(evidence["sections"], list):
            refs: set[str] = set()
            for section_id in evidence["sections"]:
                _text(section_id, where + ".sections")
                _unique(section_id, refs, where + ".sections")
                if section_id not in matching_sections:
                    _fail(where + ".sections", "cross-donor or cross-context section reference")
        else:
            _count(evidence["sections"], where + ".sections")
        _count(evidence["edges"], where + ".edges")
        if evidence["loo"] is not None:
            loo = _object(evidence["loo"], where + ".loo")
            _required(loo, ("score", "p"), where + ".loo")
            _number(loo["score"], where + ".loo.score")
            _number(loo["p"], where + ".loo.p", probability=True)


def _json(value: Any, *, pretty: bool = False) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False,
                      indent=2 if pretty else None,
                      separators=None if pretty else (",", ":"))


def _script_json(value: Any) -> str:
    # Safe even if this serialized text is subsequently embedded in HTML.
    return (_json(value).replace("&", "\\u0026").replace("<", "\\u003c")
            .replace(">", "\\u003e").replace("\u2028", "\\u2028")
            .replace("\u2029", "\\u2029"))


def _script_value(value: Any) -> str:
    # JSON.parse preserves ordinary own properties named '__proto__'; a JS
    # object literal would treat that name as a prototype setter instead.
    return "JSON.parse(" + _script_json(_json(value)) + ")"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _image_source(section: dict[str, Any]) -> Path | None:
    if section.get("image") is None:
        return None
    path = Path(section["image"]["path"]).expanduser().resolve()
    if not path.is_file():
        _fail(f"sections[{section['id']}].image.path", "image file does not exist")
    with path.open("rb") as handle:
        header = handle.read(16)
    suffix = path.suffix.lower()
    valid = ((suffix == ".png" and header.startswith(b"\x89PNG\r\n\x1a\n")) or
             (suffix in (".jpg", ".jpeg") and header.startswith(b"\xff\xd8\xff")) or
             (suffix == ".gif" and header[:6] in (b"GIF87a", b"GIF89a")) or
             (suffix == ".webp" and header[:4] == b"RIFF" and header[8:12] == b"WEBP"))
    if not valid:
        _fail(f"sections[{section['id']}].image.path", "expected a local PNG, JPEG, GIF, or WebP image")
    return path


def _tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    def cell(value: Any) -> Any:
        if value is None:
            return "null"
        if isinstance(value, (dict, list, bool)):
            return _json(value)
        # A quoted TSV field alone does not prevent spreadsheet formula execution.
        if isinstance(value, str) and value[:1] in ("=", "+", "-", "@"):
            return "'" + value
        return value
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: cell(row.get(field)) for field in fields} for row in rows)


def _columns(rows: list[dict[str, Any]], primary: list[str]) -> list[str]:
    return primary + sorted({key for row in rows for key in row} - set(primary))


def export_report(bundle: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Write a new/empty directory and return its checksum manifest.

    Relative source and image paths supplied through this Python API resolve
    against the current working directory. The CLI resolves them against the
    input JSON's directory. Existing nonempty output directories are rejected.
    """
    validate_bundle(bundle)
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise FileExistsError(f"output directory must be new or empty: {destination}")
    for name in ("index.html", "app.js", "style.css"):
        if not (_ASSETS / name).is_file():
            raise FileNotFoundError(f"bundled report asset is missing: {_ASSETS / name}")
    image_paths = [_image_source(section) for section in bundle["sections"]]
    source_manifest = []
    for source in bundle["metadata"]["sources"]:
        path = Path(source["path"]).expanduser()
        actual = _sha256(path) if path.is_file() else None
        if actual is not None and source["sha256"] is not None and actual != source["sha256"].lower():
            _fail("metadata.sources", f"SHA-256 mismatch for {source['path']!r}")
        source_manifest.append({**source, "verified_sha256": actual,
                                "verification": "verified" if actual is not None else "unavailable"})

    # Build detached presentation dictionaries; never mutate caller-owned data.
    section_summaries = []
    section_payloads = []
    image_copies: list[tuple[Path, str]] = []
    for i, (section, image_path) in enumerate(zip(bundle["sections"], image_paths)):
        file = f"data/section_{i:04d}.js"
        summary = {key: value for key, value in section.items()
                   if key not in ("spots", "edges", "expression", "image")}
        summary.update(n_spots=len(section["spots"]), n_edges=len(section["edges"]),
                       expression_genes=list(section["expression"]), file=file)
        payload = dict(section)
        if image_path is not None:
            image_file = f"assets/section_{i:04d}{image_path.suffix.lower()}"
            output_image = {**section["image"], "path": image_file}
            summary["image"] = output_image
            payload["image"] = output_image
            image_copies.append((image_path, image_file))
        section_summaries.append(summary)
        section_payloads.append((file, payload))
    report = {**bundle, "sections": section_summaries}
    index = (_ASSETS / "index.html").read_text(encoding="utf-8")
    index = index.replace("__LUCENT_TITLE__", html.escape(bundle["metadata"]["title"], quote=True))
    index = index.replace("__LUCENT_SECTION_SCRIPTS__", "")

    destination.mkdir(parents=True, exist_ok=True)
    for subdir in ("assets", "data", "downloads"):
        (destination / subdir).mkdir()
    (destination / "index.html").write_text(index, encoding="utf-8")
    for name in ("app.js", "style.css"):
        shutil.copyfile(_ASSETS / name, destination / "assets" / name)
    for source, target in image_copies:
        shutil.copyfile(source, destination / target)
    (destination / "data/report.js").write_text(
        "window.LUCENT_REPORT=" + _script_value(report) + ";\n", encoding="utf-8")
    for file, payload in section_payloads:
        # defineProperty is safe for arbitrary identifiers, including '__proto__'.
        script = ("window.LUCENT_SECTIONS=window.LUCENT_SECTIONS||Object.create(null);\n"
                  "Object.defineProperty(window.LUCENT_SECTIONS," + _script_json(payload["id"]) +
                  ",{value:" + _script_value(payload) + ",enumerable:true,configurable:true});\n")
        (destination / file).write_text(script, encoding="utf-8")
    (destination / "downloads/bundle.json").write_text(_json(bundle, pretty=True) + "\n", encoding="utf-8")

    unit_rows = []
    for unit in bundle["units"]:
        common = {key: value for key, value in unit.items() if key != "versions"}
        for version, result in unit["versions"].items():
            unit_rows.append({**common, "version": version, **result})
    _tsv(destination / "downloads/units.tsv", unit_rows, _columns(unit_rows, [
        "unit_id", "channel_id", "lr_pair", "context", "sender_state", "receiver_state",
        "version", "score", "p", "q", "status", "evaluable"]))
    for key, primary in (("genetic_links", ["channel_id", "locus_id", "variant_id", "gene_id", "pip", "pip_status"]),
                         ("donor_evidence", ["unit_id", "version", "donor_id", "score", "p", "status"]),
                         ("sections", ["id", "donor_id", "context", "label", "n_spots", "n_edges"])):
        rows = section_summaries if key == "sections" else bundle[key]
        _tsv(destination / "downloads" / f"{key}.tsv", rows, _columns(rows, primary))
    outputs = [{"path": path.relative_to(destination).as_posix(), "bytes": path.stat().st_size,
                "sha256": _sha256(path)} for path in sorted(destination.rglob("*")) if path.is_file()]
    manifest = {
        "schema_version": "1.0", "exporter": "lucent.report", "exporter_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(), "offline": True,
        "statistics_recomputed": False, "sources": source_manifest,
        "input_bundle_sha256": hashlib.sha256(_json(bundle).encode("utf-8")).hexdigest(),
        "input_bundle_hash_encoding": "compact UTF-8 JSON; insertion order; ensure_ascii=False",
        "counts": {key: len(bundle[key]) for key in ("units", "genetic_links", "donor_evidence", "sections")},
        "outputs": outputs,
        "manifest_note": "Output hashes exclude REPORT_MANIFEST.json itself.",
        "tsv_note": "Null is literal null; nested values are JSON. Formula-leading text is apostrophe-prefixed; bundle.json preserves original strings.",
    }
    (destination / "REPORT_MANIFEST.json").write_text(_json(manifest, pretty=True) + "\n", encoding="utf-8")
    return manifest


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _fail("input JSON", f"duplicate object key {key!r}")
        result[key] = value
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export existing LUCENT results as an offline report; no inference is run.")
    parser.add_argument("--input", required=True, type=Path, help="Report bundle JSON")
    parser.add_argument("--output", required=True, type=Path, help="New or empty output directory")
    args = parser.parse_args(argv)
    try:
        with args.input.open(encoding="utf-8-sig") as handle:
            bundle = json.load(handle, object_pairs_hook=_unique_pairs)
        validate_bundle(bundle)
        # Resolve declared filesystem inputs once, relative to the source bundle.
        parent = args.input.resolve().parent
        for source in bundle["metadata"]["sources"]:
            source["path"] = str((parent / source["path"]).resolve())
        for section in bundle["sections"]:
            if section.get("image") is not None:
                section["image"]["path"] = str((parent / section["image"]["path"]).resolve())
        manifest = export_report(bundle, args.output)
    except (ReportValidationError, OSError, json.JSONDecodeError) as exc:
        print(f"lucent.report: {exc}", file=sys.stderr)
        return 2
    print(f"Offline report: {args.output.resolve() / 'index.html'}")
    print(f"Exported {manifest['counts']['units']} units; statistics were not recomputed.")
    return 0
