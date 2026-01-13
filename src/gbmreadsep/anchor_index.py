from __future__ import annotations

import ast
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional


@dataclass(frozen=True)
class AnchorRecord:
    fqdn: str
    path: str
    start: int
    end: int
    hash8: str
    anchor_tag: str


def _iter_py_files(paths: Iterable[str | Path]) -> List[Path]:
    out: List[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_dir():
            out.extend(sorted(path.rglob("*.py")))
        elif path.suffix == ".py":
            out.append(path)
    return out


def _module_name_for_path(path: Path) -> str:
    parts = list(path.parts)
    if "src" in parts:
        idx = parts.index("src")
        rel = Path(*parts[idx + 1 :])
    else:
        rel = path
    return ".".join(rel.with_suffix("").parts)


def _hash_span(lines: List[str], start: int, end: int) -> str:
    cleaned = [line.strip() for line in lines[start - 1 : end]]
    digest = hashlib.sha1("\n".join(cleaned).encode("utf-8")).hexdigest()
    return digest[:8]


def _anchor_tag(path: Path, start: int, end: int, hash8: str) -> str:
    return f"{path.as_posix()}:L{start}-L{end}#{hash8}"


def build_index(
    paths: Iterable[str | Path],
    out_path: str | Path,
    *,
    include_private: bool = False,
    include_methods: bool = True,
) -> Dict[str, Dict[str, object]]:
    """Build an anchor index for public symbols in the given paths."""
    records: Dict[str, Dict[str, object]] = {}

    for py_file in _iter_py_files(paths):
        source = py_file.read_text(encoding="utf-8")
        lines = source.splitlines()
        tree = ast.parse(source)
        module = _module_name_for_path(py_file)

        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                name = node.name
                if not include_private and name.startswith("_"):
                    continue
                start, end = node.lineno, node.end_lineno
                hash8 = _hash_span(lines, start, end)
                fqdn = f"{module}.{name}"
                records[fqdn] = {
                    "path": py_file.as_posix(),
                    "start": start,
                    "end": end,
                    "hash8": hash8,
                    "anchor_tag": _anchor_tag(py_file, start, end, hash8),
                }

            if isinstance(node, ast.ClassDef):
                cname = node.name
                if not include_private and cname.startswith("_"):
                    continue
                start, end = node.lineno, node.end_lineno
                hash8 = _hash_span(lines, start, end)
                fqdn = f"{module}.{cname}"
                records[fqdn] = {
                    "path": py_file.as_posix(),
                    "start": start,
                    "end": end,
                    "hash8": hash8,
                    "anchor_tag": _anchor_tag(py_file, start, end, hash8),
                }

                if include_methods:
                    for item in node.body:
                        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                            mname = item.name
                            if not include_private and mname.startswith("_"):
                                continue
                            start, end = item.lineno, item.end_lineno
                            hash8 = _hash_span(lines, start, end)
                            fqdn = f"{module}.{cname}.{mname}"
                            records[fqdn] = {
                                "path": py_file.as_posix(),
                                "start": start,
                                "end": end,
                                "hash8": hash8,
                                "anchor_tag": _anchor_tag(py_file, start, end, hash8),
                            }

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(dict(sorted(records.items())), indent=2), encoding="utf-8")
    return records


def verify_index(index_path: str | Path) -> Dict[str, object]:
    """Verify an anchor index against current source.

    Returns a report; raises SystemExit(1) if mismatches are found.
    """
    index_path = Path(index_path)
    data = json.loads(index_path.read_text(encoding="utf-8"))

    mismatches: List[str] = []
    for fqdn, rec in data.items():
        path = Path(rec["path"])
        if not path.exists():
            mismatches.append(f"Missing file: {path}")
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        start, end = int(rec["start"]), int(rec["end"])
        hash8 = _hash_span(lines, start, end)
        if hash8 != rec["hash8"]:
            mismatches.append(f"Hash mismatch for {fqdn}: {rec['hash8']} != {hash8}")

    report = {
        "index": str(index_path),
        "symbols": len(data),
        "mismatches": mismatches,
    }

    if mismatches:
        raise SystemExit(json.dumps(report, indent=2))
    return report
