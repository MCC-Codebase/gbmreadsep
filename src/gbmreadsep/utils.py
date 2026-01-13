from __future__ import annotations

import gzip
import json
import logging
import math
import os
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, TextIO, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def logit(p: float) -> float:
    p = clamp(p, 1e-12, 1 - 1e-12)
    return math.log(p / (1 - p))


def sigmoid(x: float) -> float:
    # numerically stable sigmoid
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def phred_to_error_prob(q: int) -> float:
    # Guard against negative values (can occur if qualities are missing).
    if q <= 0:
        return 1.0
    return 10 ** (-q / 10)


def ensure_outdir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def open_textmaybe_gzip(path: str | Path, mode: str = "rt") -> TextIO:
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode)  # type: ignore[return-value]
    return open(p, mode)


def write_json(path: str | Path, obj: Any) -> None:
    with open(path, "wt", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def dataclass_to_jsonable(dc: Any) -> Mapping[str, Any]:
    return asdict(dc)


def chunked(iterable: Iterable[T], n: int) -> Iterable[list[T]]:
    chunk: list[T] = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) >= n:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
