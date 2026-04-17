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

    with pytest.raises(rcsb.RCSBFetchError) as error_info:
        rcsb._download_bytes(
            "https://example.org/missing",
            retries=5,
            retry_delay=0.0,
        )

    assert attempts["count"] == 1
    assert error_info.value.status_code == 404
    assert error_info.value.permanent is True


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


def test_parse_cluster_representative_member_entry_ids_normalizes_and_dedupes():
    cluster_text = "\n".join(
        [
            "4jpn_1 4jpp_1 4JPP_2 BADTOKEN af_test_1",
            "1abc_1 1abd_2 1ABC_3 1ABE_1",
            "not_a_pdb",
        ]
    )

    parsed = rcsb.parse_cluster_representative_member_entry_ids(cluster_text)

    assert parsed == {
        "4JPN": ["4JPN", "4JPP"],
        "1ABC": ["1ABC", "1ABD", "1ABE"],
    }


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


def test_fetch_cluster_representative_member_entry_ids_respects_cap(monkeypatch):
    cluster_text = "\n".join(
        [
            "1abc_1 1abd_1",
            "2def_1 2deg_1",
            "3ghi_1 3ghj_1",
        ]
    )
    monkeypatch.setattr(
        rcsb,
        "_download_bytes",
        lambda url, timeout=60.0, retries=0, retry_delay=1.0: cluster_text.encode("utf-8"),
    )

    parsed = rcsb.fetch_cluster_representative_member_entry_ids(30, max_structures=2)
    assert parsed == {
        "1ABC": ["1ABC", "1ABD"],
        "2DEF": ["2DEF", "2DEG"],
    }


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


def test_download_structure_cif_retries_on_corruption(monkeypatch, tmp_path: Path):
    attempts = {"count": 0}
    good_payload = b"data_TEST\n_cell.length_a 405.92\n"
    corrupt_payload = b"data_TEST\n_cell.length_a 405.\x1792\n"

    def _fake_download_bytes(url, timeout=60.0, retries=0, retry_delay=1.0):
        attempts["count"] += 1
        if attempts["count"] < 3:
            return corrupt_payload
        return good_payload

    monkeypatch.setattr(rcsb, "_download_bytes", _fake_download_bytes)
    monkeypatch.setattr(rcsb.time, "sleep", lambda seconds: None)

    output_path = rcsb.download_structure_cif(
        "1abc",
        tmp_path,
        retry_delay=0.0,
        corruption_retries=4,
    )
    assert output_path.read_bytes() == good_payload
    assert attempts["count"] == 3


def test_download_structure_cif_raises_on_persistent_corruption(monkeypatch, tmp_path: Path):
    attempts = {"count": 0}
    corrupt_payload = b"data_TEST\n_cell.length_a 405.\x1792\n"

    def _fake_download_bytes(url, timeout=60.0, retries=0, retry_delay=1.0):
        attempts["count"] += 1
        return corrupt_payload

    monkeypatch.setattr(rcsb, "_download_bytes", _fake_download_bytes)
    monkeypatch.setattr(rcsb.time, "sleep", lambda seconds: None)

    with pytest.raises(rcsb.RCSBFetchError) as error_info:
        rcsb.download_structure_cif(
            "1abc",
            tmp_path,
            retry_delay=0.0,
            corruption_retries=4,
        )

    assert attempts["count"] == 5
    assert error_info.value.permanent is False
    assert not (tmp_path / "1ABC.cif").exists()


def test_find_cif_corruption_reason_detects_control_bytes_and_empty():
    assert rcsb._find_cif_corruption_reason(b"") == "empty payload"
    assert rcsb._find_cif_corruption_reason(b"data_OK\n405.92\n") is None
    assert rcsb._find_cif_corruption_reason(b"405.\x1792") == "contains control byte 0x17"
    # tab / newline / carriage return are fine
    assert rcsb._find_cif_corruption_reason(b"a\tb\nc\rd") is None
