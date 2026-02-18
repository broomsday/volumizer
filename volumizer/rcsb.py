"""
Helpers for downloading structures and sequence clusters from RCSB services.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


RCSB_CIF_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_CLUSTER_URL_TEMPLATE = (
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{identity}.txt"
)
RCSB_ENTRY_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

VALID_CLUSTER_IDENTITIES = frozenset([30, 40, 50, 70, 90, 95, 100])
TERMINAL_HTTP_STATUS_CODES = frozenset([400, 401, 403, 404])

EXPERIMENTAL_METHOD_FILTERS = {
    "xray": frozenset(["X-RAY DIFFRACTION"]),
    "em": frozenset(["ELECTRON MICROSCOPY"]),
    "nmr": frozenset(["SOLUTION NMR", "SOLID-STATE NMR"]),
    "neutron": frozenset(["NEUTRON DIFFRACTION"]),
}

VALID_CLUSTER_METHOD_FILTERS = frozenset(EXPERIMENTAL_METHOD_FILTERS.keys())
DEFAULT_CLUSTER_METHOD_FILTERS = ("xray", "em")


class RCSBFetchError(RuntimeError):
    """
    Structured fetch failure for RCSB network calls.
    """

    def __init__(
        self,
        message: str,
        *,
        url: str,
        status_code: int | None = None,
        reason: str | None = None,
        permanent: bool = False,
    ):
        super().__init__(message)
        self.url = url
        self.status_code = status_code
        self.reason = reason
        self.permanent = bool(permanent)


def normalize_pdb_id(pdb_id: str) -> str:
    """
    Normalize and validate a 4-character PDB ID.
    """
    normalized = pdb_id.strip().upper()
    if len(normalized) != 4 or not normalized.isalnum():
        raise ValueError(f"Invalid PDB ID: {pdb_id!r}")
    return normalized


def is_pdb_entry_id(entry_id: str) -> bool:
    """
    Return True for 4-character alphanumeric PDB entry identifiers.
    """
    return len(entry_id) == 4 and entry_id.isalnum()


def normalize_method_filter_name(method_filter: str) -> str:
    """
    Normalize CLI method aliases to canonical filter names.
    """
    normalized = method_filter.strip().lower().replace("_", "-")
    alias_map = {
        "x-ray": "xray",
        "xray": "xray",
        "em": "em",
        "cryo-em": "em",
        "cryoem": "em",
        "electron-microscopy": "em",
        "nmr": "nmr",
        "neutron": "neutron",
    }
    canonical = alias_map.get(normalized, normalized)
    if canonical not in VALID_CLUSTER_METHOD_FILTERS:
        allowed = ", ".join(sorted(VALID_CLUSTER_METHOD_FILTERS))
        raise ValueError(
            f"Unsupported cluster method filter {method_filter!r}. Allowed: {allowed}."
        )
    return canonical


def _download_bytes(
    url: str,
    timeout: float = 60.0,
    retries: int = 0,
    retry_delay: float = 1.0,
) -> bytes:
    """
    Download raw bytes from a URL with simple exponential backoff retries.
    """
    max_retries = max(0, int(retries))
    base_delay = max(0.0, float(retry_delay))

    for attempt in range(max_retries + 1):
        try:
            with urlopen(url, timeout=timeout) as response:
                return response.read()
        except HTTPError as error:
            status_code = int(error.code)
            is_terminal = status_code in TERMINAL_HTTP_STATUS_CODES
            if attempt >= max_retries or is_terminal:
                raise RCSBFetchError(
                    f"HTTP error while fetching {url}: {status_code}",
                    url=url,
                    status_code=status_code,
                    reason=error.reason,
                    permanent=is_terminal,
                ) from error
        except URLError as error:
            if attempt >= max_retries:
                raise RCSBFetchError(
                    f"Network error while fetching {url}: {error.reason}",
                    url=url,
                    reason=str(error.reason),
                    permanent=False,
                ) from error

        if base_delay > 0:
            time.sleep(base_delay * (2**attempt))

    raise RCSBFetchError(
        f"Failed to fetch {url}",
        url=url,
        permanent=False,
    )


def download_structure_cif(
    pdb_id: str,
    output_dir: Path,
    overwrite: bool = False,
    timeout: float = 60.0,
    retries: int = 0,
    retry_delay: float = 1.0,
) -> Path:
    """
    Download a structure CIF by PDB ID and return local file path.
    """
    normalized_id = normalize_pdb_id(pdb_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    out_path = output_dir / f"{normalized_id}.cif"
    if out_path.is_file() and not overwrite:
        return out_path

    url = RCSB_CIF_URL_TEMPLATE.format(pdb_id=normalized_id)
    data = _download_bytes(
        url,
        timeout=timeout,
        retries=retries,
        retry_delay=retry_delay,
    )
    out_path.write_bytes(data)
    return out_path


def fetch_entry_metadata(
    pdb_id: str,
    timeout: float = 60.0,
    retries: int = 0,
    retry_delay: float = 1.0,
) -> dict:
    """
    Fetch RCSB core entry JSON metadata for one PDB ID.
    """
    normalized_id = normalize_pdb_id(pdb_id)
    url = RCSB_ENTRY_URL_TEMPLATE.format(pdb_id=normalized_id)
    payload = _download_bytes(
        url,
        timeout=timeout,
        retries=retries,
        retry_delay=retry_delay,
    )
    return json.loads(payload.decode("utf-8"))


def extract_entry_experimental_methods(entry_metadata: dict) -> list[str]:
    """
    Extract normalized experimental method labels from entry metadata.
    """
    methods: list[str] = []
    seen: set[str] = set()

    for exptl_item in entry_metadata.get("exptl", []):
        if not isinstance(exptl_item, dict):
            continue
        method_value = exptl_item.get("method")
        if not isinstance(method_value, str):
            continue
        normalized_method = method_value.strip().upper()
        if len(normalized_method) == 0 or normalized_method in seen:
            continue
        methods.append(normalized_method)
        seen.add(normalized_method)

    return methods


def extract_entry_best_resolution(entry_metadata: dict) -> float | None:
    """
    Return the best (lowest) combined resolution, if present.
    """
    entry_info = entry_metadata.get("rcsb_entry_info")
    if not isinstance(entry_info, dict):
        return None

    resolution_values = entry_info.get("resolution_combined")
    if not isinstance(resolution_values, list):
        return None

    numeric_values = []
    for value in resolution_values:
        try:
            numeric_values.append(float(value))
        except (TypeError, ValueError):
            continue

    if len(numeric_values) == 0:
        return None
    return min(numeric_values)


def _expand_allowed_methods(allowed_method_filters: list[str] | tuple[str, ...]) -> set[str]:
    methods: set[str] = set()
    for method_filter in allowed_method_filters:
        canonical_filter = normalize_method_filter_name(method_filter)
        methods.update(EXPERIMENTAL_METHOD_FILTERS[canonical_filter])
    return methods


def entry_passes_filters(
    entry_metadata: dict,
    allowed_method_filters: list[str] | tuple[str, ...] | None = None,
    max_resolution: float | None = None,
) -> tuple[bool, str | None]:
    """
    Apply experimental-method and resolution filters to one entry.

    Returns `(passes, failure_reason)` where `failure_reason` is one of:
    - `experimental_method`
    - `resolution`
    - `missing_resolution`
    """
    if allowed_method_filters is not None:
        entry_methods = set(extract_entry_experimental_methods(entry_metadata))
        allowed_methods = _expand_allowed_methods(allowed_method_filters)
        if len(entry_methods.intersection(allowed_methods)) == 0:
            return False, "experimental_method"

    if max_resolution is not None:
        best_resolution = extract_entry_best_resolution(entry_metadata)
        if best_resolution is None:
            return False, "missing_resolution"
        if best_resolution > max_resolution:
            return False, "resolution"

    return True, None


def parse_cluster_representative_entry_ids(
    cluster_text: str,
    include_non_pdb: bool = False,
) -> list[str]:
    """
    Parse sequence cluster text and return de-duplicated representative entry IDs.

    Cluster lines contain polymer entity IDs, e.g. `1ABC_1 2XYZ_2 ...`.
    The first token on each line is the cluster representative.
    """
    representatives = []
    seen = set()

    for raw_line in cluster_text.splitlines():
        line = raw_line.strip()
        if len(line) == 0:
            continue
        tokens = line.split()
        if len(tokens) == 0:
            continue

        representative_entity = tokens[0]
        entry_id = representative_entity.split("_", 1)[0].upper()

        if not include_non_pdb and not is_pdb_entry_id(entry_id):
            continue
        if entry_id in seen:
            continue

        representatives.append(entry_id)
        seen.add(entry_id)

    return representatives


def fetch_cluster_representative_entry_ids(
    identity: int,
    max_structures: int | None = None,
    timeout: float = 60.0,
    include_non_pdb: bool = False,
    retries: int = 0,
    retry_delay: float = 1.0,
) -> list[str]:
    """
    Fetch and parse representative entry IDs for a sequence-identity cluster file.
    """
    if identity not in VALID_CLUSTER_IDENTITIES:
        allowed = ", ".join([str(value) for value in sorted(VALID_CLUSTER_IDENTITIES)])
        raise ValueError(
            f"Unsupported identity threshold {identity}. Allowed: {allowed}."
        )

    url = RCSB_CLUSTER_URL_TEMPLATE.format(identity=identity)
    cluster_text = _download_bytes(
        url,
        timeout=timeout,
        retries=retries,
        retry_delay=retry_delay,
    ).decode("utf-8")
    entry_ids = parse_cluster_representative_entry_ids(
        cluster_text,
        include_non_pdb=include_non_pdb,
    )

    if max_structures is not None and max_structures > 0:
        return entry_ids[:max_structures]

    return entry_ids
