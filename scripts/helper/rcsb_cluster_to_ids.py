"""
Open an RCSB cluster file and generate a text file with one ID per line,
where each ID is the first ID in the cluster.
"""


from pathlib import Path

import typer


def main(cluster_file: Path, output_file: Path) -> None:
    """
    Open an RCSB cluster file and generate a text file with one ID per line,
    where each ID is the first ID in the cluster.
    """

    with open(cluster_file, mode="r", encoding="utf-8") as file_in:
        pdbs = {f"{line.split()[0][:4]}\n" for line in file_in.readlines()}

    with open(output_file, mode="w", encoding="utf-8") as file_out:
        file_out.writelines(pdbs)


if "__main__" in __name__:
    typer.run(main)