"""
Helpers for downloading structures and sequence clusters from RCSB services.
"""

from __future__ import annotations

from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


RCSB_CIF_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_CLUSTER_URL_TEMPLATE = (
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{identity}.txt"
)

VALID_CLUSTER_IDENTITIES = frozenset([30, 40, 50, 70, 90, 95, 100])


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


def _download_bytes(url: str, timeout: float = 60.0) -> bytes:
    """
    Download raw bytes from a URL.
    """
    try:
        with urlopen(url, timeout=timeout) as response:
            return response.read()
    except HTTPError as error:
        raise RuntimeError(f"HTTP error while fetching {url}: {error.code}") from error
    except URLError as error:
        raise RuntimeError(f"Network error while fetching {url}: {error.reason}") from error


def download_structure_cif(
    pdb_id: str,
    output_dir: Path,
    overwrite: bool = False,
    timeout: float = 60.0,
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
    data = _download_bytes(url, timeout=timeout)
    out_path.write_bytes(data)
    return out_path


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
    cluster_text = _download_bytes(url, timeout=timeout).decode("utf-8")
    entry_ids = parse_cluster_representative_entry_ids(
        cluster_text,
        include_non_pdb=include_non_pdb,
    )

    if max_structures is not None and max_structures > 0:
        return entry_ids[:max_structures]

    return entry_ids
