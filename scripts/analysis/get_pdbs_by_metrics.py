"""
Scan over annotated DFs for pores, pockets, and/or cavities matching given
metric constraints.

Example:
    python scripts/analysis/get_pdbs_by_metrics.py data/rcsb_cluster/bc-90.txt pore_analysis_rcsb_3.0.txt --find-pores --min-dimension-one=30.0 --min-dimension-two=15.0 --min-dimension-three=15.0
"""


from pathlib import Path
import warnings

import typer

from pore import analysis
from pore.cli import guess_analysis_input_type


def main(
    analysis_input: Path = typer.Argument(..., help=""),
    analysis_output: Path = typer.Argument(..., help=""),
    find_pores: bool = typer.Option(False, help=""),
    find_pockets: bool = typer.Option(False, help=""),
    find_cavities: bool = typer.Option(False, help=""),
    min_volume: float = typer.Option(0.0, help=""),
    max_volume: float = typer.Option(None, help=""),
    min_dimension_one: float = typer.Option(0.0, help=""),
    max_dimension_one: float = typer.Option(None, help=""),
    min_dimension_two: float = typer.Option(0.0, help=""),
    max_dimension_two: float = typer.Option(None, help=""),
    min_dimension_three: float = typer.Option(0.0, help=""),
    max_dimension_three: float = typer.Option(None, help=""),
):
    """
    Scan over annotated DFs for pores, pockets, and/or cavities matching given
    metric constraints.
    """
    if (not find_pores) and (not find_pockets) and (not find_cavities):
        warnings.warn("You have not selected any volume types to find!")

    input_type = guess_analysis_input_type(analysis_input)
    if input_type == "file":
        with open(analysis_input, mode="r", encoding="utf-8") as id_file:
            pdb_ids = [line.strip() for line in id_file.readlines()]
        annotation_paths, missing_dfs = analysis.get_annotations_by_id(pdb_ids)
        if missing_dfs > 0:
            warnings.warn(
                f"Missing {missing_dfs} dataframes",
            )
    elif input_type == "dir":
        annotation_paths = list(Path(analysis_input).glob("*.json"))
    else:
        raise RuntimeError("Input type not implemented")

    pdb_annotations = analysis.get_pdb_annotations(annotation_paths)

    metrics = {
        "pores": find_pores,
        "pockets": find_pockets,
        "cavities": find_cavities,
        "min_volume": min_volume,
        "max_volume": max_volume,
        "min_x": min_dimension_one,
        "max_x": max_dimension_one,
        "min_y": min_dimension_two,
        "max_y": max_dimension_two,
        "min_z": min_dimension_three,
        "max_z": max_dimension_three,
    }
    selected_pdb_annotations = analysis.select_annotations_by_metrics(pdb_annotations, metrics)

    annotation_names = [f"{name}\n" for name in selected_pdb_annotations.keys()]
    with open(analysis_output, mode="w", encoding="utf-8") as out_file:
        out_file.writelines(annotation_names)

    print(f"Found {len(annotation_names)} matching PDBs")


if "__main__" in __name__:
    typer.run(main)