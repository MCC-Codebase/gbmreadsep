"""Environment self-checks.

This module powers the ``gbmreadsep doctor`` CLI command.

Rationale
---------
Most GBMReadSep functionality is pure Python. However, the optional Nanopore
end-to-end workflow relies on external tools (e.g. minimap2, samtools, Docker).
For biologists with minimal computational background, a single command that
pinpoints missing dependencies is often the most helpful UX.
"""

from __future__ import annotations

import logging
import platform
import shutil
from dataclasses import dataclass
from typing import Dict, Optional

from .external import run_command

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CheckResult:
    name: str
    ok: bool
    detail: str
    howto: Optional[str] = None


def _which(name: str) -> Optional[str]:
    return shutil.which(name)


def check_python() -> CheckResult:
    v = platform.python_version()
    return CheckResult(name="python", ok=True, detail=f"Python {v}")


def check_executable(name: str, *, howto: Optional[str] = None) -> CheckResult:
    p = _which(name)
    if p is None:
        return CheckResult(name=name, ok=False, detail="not found in PATH", howto=howto)
    return CheckResult(name=name, ok=True, detail=p)


def check_docker() -> CheckResult:
    p = _which("docker")
    if p is None:
        return CheckResult(
            name="docker",
            ok=False,
            detail="not found in PATH",
            howto=(
                "Install Docker Desktop (Windows/Mac) or docker-ce (Linux) and ensure the 'docker' command works.\n"
                "Then re-run: gbmreadsep doctor"
            ),
        )
    try:
        cp = run_command(["docker", "version"], check=True, capture=True, text=True)
        # Keep short detail; docker version output is verbose.
        detail = "docker OK"
        if isinstance(cp.stdout, str) and cp.stdout.strip():
            detail = "docker OK (client/server reachable)"
        return CheckResult(name="docker", ok=True, detail=detail)
    except Exception as e:
        return CheckResult(
            name="docker",
            ok=False,
            detail=f"docker present but not usable: {e}",
            howto=(
                "Ensure the Docker daemon is running and you have permission to access it.\n"
                "On Linux, you may need to add your user to the 'docker' group or use sudo."
            ),
        )


def collect_checks() -> Dict[str, CheckResult]:
    """Run all checks and return a mapping name->result."""
    checks: Dict[str, CheckResult] = {}

    checks["python"] = check_python()
    checks["minimap2"] = check_executable(
        "minimap2",
        howto=(
            "Ubuntu: sudo apt-get install -y minimap2\n"
            "Conda/mamba: mamba install -c bioconda minimap2"
        ),
    )
    checks["samtools"] = check_executable(
        "samtools",
        howto=(
            "Ubuntu: sudo apt-get install -y samtools\n"
            "Conda/mamba: mamba install -c bioconda samtools"
        ),
    )
    checks["bcftools"] = check_executable(
        "bcftools",
        howto=(
            "Ubuntu: sudo apt-get install -y bcftools\n"
            "Conda/mamba: mamba install -c bioconda bcftools"
        ),
    )
    checks["tabix"] = check_executable(
        "tabix",
        howto=(
            "Ubuntu: sudo apt-get install -y tabix\n"
            "Conda/mamba: mamba install -c bioconda tabix"
        ),
    )
    checks["pigz"] = check_executable(
        "pigz",
        howto=(
            "Ubuntu: sudo apt-get install -y pigz\n"
            "Conda/mamba: mamba install -c conda-forge pigz"
        ),
    )
    checks["docker"] = check_docker()

    return checks
