from pathlib import Path

import pytest

from gbmreadsep.anchor_index import build_index, verify_index


def test_anchor_index_build_and_verify(tmp_path: Path) -> None:
    src = tmp_path / "sample.py"
    src.write_text(
        """
class Demo:\n    def method(self):\n        return 1\n\n\ndef public_fn():\n    return 2\n""".lstrip(),
        encoding="utf-8",
    )

    index_path = tmp_path / "index.json"
    build_index([src], index_path)

    report = verify_index(index_path)
    assert report["symbols"] >= 2

    # Modify source to cause mismatch
    src.write_text(
        """
class Demo:\n    def method(self):\n        return 3\n\n\ndef public_fn():\n    return 2\n""".lstrip(),
        encoding="utf-8",
    )

    with pytest.raises(SystemExit):
        verify_index(index_path)
