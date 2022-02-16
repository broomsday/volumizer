"""
Scan over annotated DFs for pores, pockets, and/or cavities matching given
metric constraints.
"""


from pathlib import Path
import warnings

import typer

from pore import analysis
from pore.cli import guess_analysis_input_type


def main(
    analysis_input: str = typer.Argument(..., help=""),
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
    selected_pdb_annotations = analysis.select_annotations_by_type(
        pdb_annotations, find_pores, find_pockets, find_cavities
    )

    metrics = {
        "min_volume": min_volume,
        "max_volume": max_volume,
        "min_dimension_one": min_dimension_one,
        "max_dimension_one": max_dimension_one,
    }
    selected_pdb_annotations = analysis.select_annotations_by_metrics(selected_pdb_annotations, metrics)

    print(selected_pdb_annotations[list(pdb_annotations.keys())[0]])
    print(len(selected_pdb_annotations))


if "__main__" in __name__:
    typer.run(main)