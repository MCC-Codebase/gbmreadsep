"""Helpers for running external commands (minimap2/samtools/docker/etc).

Design goals
------------
- Fail fast with actionable error messages.
- Capture stderr/stdout for debugging.
- Keep the public surface small; treat this as an internal utility module.

This package intentionally *does not* vendor heavyweight bioinformatics tooling.
Instead, it provides robust wrappers around well-established external tools.
"""

from __future__ import annotations

import logging
import os
import shlex
import subprocess
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)


class ExternalCommandError(RuntimeError):
    """Raised when an external command fails."""

    def __init__(
        self,
        message: str,
        *,
        cmd: Sequence[str],
        returncode: int,
        stdout: Optional[str] = None,
        stderr: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.cmd = list(cmd)
        self.returncode = int(returncode)
        self.stdout = stdout
        self.stderr = stderr


def cmd_to_str(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


def ensure_executable_in_path(exe: str, *, hint: Optional[str] = None) -> None:
    """Ensure an executable exists in PATH.

    Parameters
    ----------
    exe:
        Name of the executable to find.
    hint:
        Optional message shown if the executable is missing.
    """
    from shutil import which

    if which(exe) is None:
        msg = f"Required executable '{exe}' was not found in your PATH."
        if hint:
            msg += "\n\n" + hint
        raise FileNotFoundError(msg)


def run_command(
    cmd: Sequence[str],
    *,
    cwd: Optional[str | Path] = None,
    env: Optional[Mapping[str, str]] = None,
    check: bool = True,
    capture: bool = True,
    text: bool = True,
) -> subprocess.CompletedProcess:
    """Run a command and return the CompletedProcess.

    If ``check`` is True, raise ``ExternalCommandError`` on non-zero exit.

    Notes
    -----
    - We default to capturing stdout+stderr to improve error messages.
    - For large streaming pipelines (e.g. minimap2 | samtools sort), use
      :func:`run_pipe` instead.
    """
    if cwd is not None:
        cwd = str(Path(cwd))

    env_merged: Optional[Dict[str, str]]
    if env is None:
        env_merged = None
    else:
        env_merged = dict(os.environ)
        env_merged.update({str(k): str(v) for k, v in env.items()})

    logger.debug("Running command: %s", cmd_to_str(cmd))

    cp = subprocess.run(
        list(map(str, cmd)),
        cwd=cwd,
        env=env_merged,
        check=False,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        text=text,
    )

    if check and cp.returncode != 0:
        raise ExternalCommandError(
            message=textwrap.dedent(
                f"""
                External command failed (exit code {cp.returncode}).

                Command:
                  {cmd_to_str(cmd)}

                STDERR (tail):
                  {(_tail(cp.stderr) if isinstance(cp.stderr, str) else str(cp.stderr))}
                """
            ).strip(),
            cmd=cmd,
            returncode=cp.returncode,
            stdout=cp.stdout if isinstance(cp.stdout, str) else None,
            stderr=cp.stderr if isinstance(cp.stderr, str) else None,
        )

    return cp


def _tail(s: Optional[str], n: int = 3000) -> str:
    if not s:
        return "(empty)"
    s = str(s)
    if len(s) <= n:
        return s
    return "..." + s[-n:]


def run_pipe(
    producer: Sequence[str],
    consumer: Sequence[str],
    *,
    cwd: Optional[str | Path] = None,
    env: Optional[Mapping[str, str]] = None,
    check: bool = True,
) -> Tuple[subprocess.CompletedProcess, subprocess.CompletedProcess]:
    """Run a two-process pipe: ``producer | consumer``.

    Returns
    -------
    (cp1, cp2)
        CompletedProcess-like objects (stdout/stderr captured).

    Raises
    ------
    ExternalCommandError
        If ``check`` is True and any command fails.
    """
    if cwd is not None:
        cwd = str(Path(cwd))

    env_merged: Optional[Dict[str, str]]
    if env is None:
        env_merged = None
    else:
        env_merged = dict(os.environ)
        env_merged.update({str(k): str(v) for k, v in env.items()})

    logger.debug("Running pipe: %s | %s", cmd_to_str(producer), cmd_to_str(consumer))

    p1 = subprocess.Popen(
        list(map(str, producer)),
        cwd=cwd,
        env=env_merged,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,  # bytes to avoid encoding overhead
    )
    p2 = subprocess.Popen(
        list(map(str, consumer)),
        cwd=cwd,
        env=env_merged,
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )

    assert p1.stdout is not None
    p1.stdout.close()  # allow p1 to receive SIGPIPE if p2 exits

    out2, err2 = p2.communicate()
    out1, err1 = p1.communicate()

    # Build CompletedProcess wrappers
    cp1 = subprocess.CompletedProcess(args=list(producer), returncode=p1.returncode, stdout=out1, stderr=err1)
    cp2 = subprocess.CompletedProcess(args=list(consumer), returncode=p2.returncode, stdout=out2, stderr=err2)

    if check and (cp1.returncode != 0 or cp2.returncode != 0):
        # Prefer the consumer stderr, but include both tails.
        msg = textwrap.dedent(
            f"""
            External pipe failed.

            Producer:
              {cmd_to_str(producer)}
              exit={cp1.returncode}
              stderr tail: {_tail_bytes(cp1.stderr)}

            Consumer:
              {cmd_to_str(consumer)}
              exit={cp2.returncode}
              stderr tail: {_tail_bytes(cp2.stderr)}
            """
        ).strip()
        raise ExternalCommandError(
            msg,
            cmd=consumer,
            returncode=cp2.returncode if cp2.returncode != 0 else cp1.returncode,
            stdout=None,
            stderr=None,
        )

    return cp1, cp2


def _tail_bytes(b: Optional[bytes], n: int = 3000) -> str:
    if not b:
        return "(empty)"
    try:
        s = b.decode("utf-8", errors="replace")
    except Exception:
        s = repr(b)
    return _tail(s, n=n)


@dataclass(frozen=True)
class DockerMount:
    host_dir: Path
    container_dir: str
    read_only: bool = True


def build_docker_mounts(
    host_paths: Iterable[str | Path],
    *,
    container_root: str = "/mnt/gbmreadsep",
) -> Tuple[List[DockerMount], Dict[Path, str]]:
    """Create a minimal set of directory mounts for a set of host file paths.

    Parameters
    ----------
    host_paths:
        File or directory paths on the host.

    Returns
    -------
    mounts:
        Mount specifications (unique host directories only).
    mapping:
        Mapping from host directory -> container directory.

    Notes
    -----
    We mount *directories*, not individual files, because Docker cannot mount files
    across platforms consistently.
    """
    uniq_dirs: List[Path] = []
    for p in host_paths:
        hp = Path(p).expanduser().resolve()
        d = hp if hp.is_dir() else hp.parent
        if d not in uniq_dirs:
            uniq_dirs.append(d)

    mapping: Dict[Path, str] = {}
    mounts: List[DockerMount] = []

    for i, d in enumerate(uniq_dirs):
        cdir = f"{container_root}/in{i}"
        mapping[d] = cdir
        mounts.append(DockerMount(host_dir=d, container_dir=cdir, read_only=True))

    return mounts, mapping


def docker_run(
    *,
    image: str,
    argv: Sequence[str],
    mounts: Sequence[DockerMount],
    workdir: Optional[str] = None,
    env: Optional[Mapping[str, str]] = None,
    extra_docker_args: Optional[Sequence[str]] = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run a docker container and execute argv inside.

    Parameters
    ----------
    image:
        Docker image, e.g. ``hkubal/clairs:latest``.
    argv:
        Command executed inside the container.
    mounts:
        List of mounts to add (directories).
    workdir:
        Optional working directory inside the container.
    env:
        Environment variables for the container.
    extra_docker_args:
        Extra args inserted right after ``docker run``.

    Returns
    -------
    CompletedProcess
    """
    ensure_executable_in_path(
        "docker",
        hint="Install Docker Desktop (Windows/Mac) or docker-ce (Linux), then ensure 'docker' works from your shell.",
    )

    cmd: List[str] = ["docker", "run", "--rm"]

    if extra_docker_args:
        cmd.extend(list(map(str, extra_docker_args)))

    for m in mounts:
        ro = ":ro" if m.read_only else ""
        cmd.extend(["-v", f"{str(m.host_dir)}:{m.container_dir}{ro}"])

    if env:
        for k, v in env.items():
            cmd.extend(["-e", f"{k}={v}"])

    if workdir:
        cmd.extend(["-w", workdir])

    cmd.append(image)
    cmd.extend(list(map(str, argv)))

    return run_command(cmd, check=check, capture=True, text=True)
