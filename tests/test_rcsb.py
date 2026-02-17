from pathlib import Path
from urllib.error import HTTPError, URLError

import pytest

from volumizer import rcsb


def test_normalize_pdb_id_valid_and_invalid():
    assert rcsb.normalize_pdb_id("4jpn") == "4JPN"

    with pytest.raises(ValueError):
        rcsb.normalize_pdb_id("bad")
    with pytest.raises(ValueError):
        rcsb.normalize_pdb_id("ABCDE")


def test_normalize_method_filter_name_aliases_and_invalid():
    assert rcsb.normalize_method_filter_name("xray") == "xray"
    assert rcsb.normalize_method_filter_name("x-ray") == "xray"
    assert rcsb.normalize_method_filter_name("cryo-em") == "em"
    assert rcsb.normalize_method_filter_name("em") == "em"

    with pytest.raises(ValueError):
        rcsb.normalize_method_filter_name("sas")


def test_entry_passes_filters_by_method_and_resolution():
    entry_metadata = {
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "rcsb_entry_info": {"resolution_combined": [2.1]},
    }

    assert rcsb.entry_passes_filters(entry_metadata, ["xray"], None) == (True, None)
    assert rcsb.entry_passes_filters(entry_metadata, ["em"], None) == (
        False,
        "experimental_method",
    )
    assert rcsb.entry_passes_filters(entry_metadata, ["xray"], 2.0) == (
        False,
        "resolution",
    )

    no_resolution_metadata = {
        "exptl": [{"method": "ELECTRON MICROSCOPY"}],
        "rcsb_entry_info": {},
    }
    assert rcsb.entry_passes_filters(no_resolution_metadata, ["em"], 3.0) == (
        False,
        "missing_resolution",
    )


def test_download_bytes_retries_then_succeeds(monkeypatch):
    attempts = {"count": 0}

    class _DummyResponse:
        def __init__(self, payload: bytes):
            self._payload = payload

        def read(self) -> bytes:
            return self._payload

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    def _fake_urlopen(url, timeout=60.0):
        attempts["count"] += 1
        if attempts["count"] < 3:
            raise URLError("temporary")
        return _DummyResponse(b"ok")

    sleep_calls = []
    monkeypatch.setattr(rcsb, "urlopen", _fake_urlopen)
    monkeypatch.setattr(rcsb.time, "sleep", lambda seconds: sleep_calls.append(seconds))

    payload = rcsb._download_bytes(
        "https://example.org/test",
        retries=2,
        retry_delay=0.5,
    )

    assert payload == b"ok"
    assert attempts["count"] == 3
    assert sleep_calls == [0.5, 1.0]


def test_download_bytes_does_not_retry_terminal_http_error(monkeypatch):
    attempts = {"count": 0}

    def _fake_urlopen(url, timeout=60.0):
        attempts["count"] += 1
        raise HTTPError(url=url, code=404, msg="not found", hdrs=None, fp=None)

    monkeypatch.setattr(rcsb, "urlopen", _fake_urlopen)

    with pytest.raises(RuntimeError):
        rcsb._download_bytes(
            "https://example.org/missing",
            retries=5,
            retry_delay=0.0,
        )

    assert attempts["count"] == 1


def test_parse_cluster_representative_entry_ids_filters_non_pdb_and_dedupes():
    cluster_text = "\n".join(
        [
            "1abc_1 9xyz_1",
            "2def_2 1abc_1",
            "AF_TEST_1 4jkl_1",
            "1ABC_3 4jkl_2",
        ]
    )

    parsed = rcsb.parse_cluster_representative_entry_ids(cluster_text)
    assert parsed == ["1ABC", "2DEF"]

    parsed_with_non_pdb = rcsb.parse_cluster_representative_entry_ids(
        cluster_text,
        include_non_pdb=True,
    )
    assert parsed_with_non_pdb == ["1ABC", "2DEF", "AF"]


def test_fetch_cluster_representative_entry_ids_respects_cap(monkeypatch):
    cluster_text = "\n".join(
        [
            "1abc_1 9xyz_1",
            "2def_2 1abc_1",
            "3ghi_1 5jkl_1",
        ]
    )
    monkeypatch.setattr(
        rcsb,
        "_download_bytes",
        lambda url, timeout=60.0, retries=0, retry_delay=1.0: cluster_text.encode("utf-8"),
    )

    parsed = rcsb.fetch_cluster_representative_entry_ids(30, max_structures=2)
    assert parsed == ["1ABC", "2DEF"]


def test_download_structure_cif_writes_downloaded_data(monkeypatch, tmp_path: Path):
    expected_content = b"data_TEST"
    monkeypatch.setattr(
        rcsb,
        "_download_bytes",
        lambda url, timeout=60.0, retries=0, retry_delay=1.0: expected_content,
    )

    output_path = rcsb.download_structure_cif("1abc", tmp_path)
    assert output_path == tmp_path / "1ABC.cif"
    assert output_path.read_bytes() == expected_content

    # Cached file should be reused when overwrite=False.
    monkeypatch.setattr(
        rcsb,
        "_download_bytes",
        lambda url, timeout=60.0, retries=0, retry_delay=1.0: (_ for _ in ()).throw(
            RuntimeError("should not download")
        ),
    )
    cached_path = rcsb.download_structure_cif("1ABC", tmp_path, overwrite=False)
    assert cached_path == output_path
