import subprocess
import sys


def test_cli_help() -> None:
    cp = subprocess.run(
        [sys.executable, "-m", "gbmreadsep", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "GBMReadSep" in cp.stdout or "gbmreadsep" in cp.stdout.lower()
