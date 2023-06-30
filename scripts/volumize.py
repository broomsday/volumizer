"""
Command-line entry-point to find pores and cavities in PDBs.
"""

from pathlib import Path

import typer
import multiprocessing

from volumizer import cli, utils, pdb, volumizer, rcsb
from volumizer.constants import VOXEL_SIZE
from volumizer.paths import PREPARED_PDB_DIR, ANNOTATED_DF_DIR, ANNOTATED_PDB_DIR


def volumize_pdb_id(pdb_id: str) -> None:
    """
    Download the given PDB ID and then volumize it.
    """
    if not utils.have_annotation(pdb_id):
        pdb_file = rcsb.download_pdb_file(pdb_id)
        if pdb_file is None:
            return None

        print(f"Working on: {pdb_id}")

        pdb_structure = pdb.load_structure(pdb_file)
        prepared_structure, preparation_remarks = volumizer.prepare_pdb_structure(pdb_structure)
        prepared_pdb = pdb.structure_to_lines(prepared_structure)
        prepared_lines = pdb.merge_pdb_lines([preparation_remarks, prepared_pdb])
        pdb.save_pdb_lines(prepared_lines, (PREPARED_PDB_DIR / pdb_file.stem).with_suffix('.pdb'))

        annotation, annotated_lines = volumizer.annotate_pdb_structure(prepared_structure)
        final_pdb_lines = pdb.merge_pdb_lines([prepared_lines, ["END"], annotated_lines])
        annotation_df = utils.make_annotation_dataframe(annotation)

        utils.save_annotation_dataframe(annotation_df, (ANNOTATED_DF_DIR / pdb_file.stem).with_suffix('.json'))
        print(f"Annotation saved to dataframe as: {(ANNOTATED_DF_DIR / pdb_file.stem).with_suffix('.json')}")
        pdb.save_pdb_lines(final_pdb_lines, (ANNOTATED_PDB_DIR / pdb_file.stem).with_suffix('.pdb'))
        print(f"Annotation saved to PDB as: {(ANNOTATED_PDB_DIR / pdb_file.stem).with_suffix('.pdb')}")
        print(f"Quick annotation output:")
        print(annotation_df)


def volumize_pdb_file(pdb_file: Path) -> None:
    """
    Operate directly on the given PDB file.
    """
    if not utils.have_annotation(pdb_file.stem):
        print(f"Working on: {pdb_file.stem}")

        pdb_structure = pdb.load_structure(pdb_file)
        prepared_structure, preparation_remarks = volumizer.prepare_pdb_structure(pdb_structure)
        prepared_pdb = pdb.structure_to_lines(prepared_structure)
        prepared_lines = pdb.merge_pdb_lines([preparation_remarks, prepared_pdb])
        pdb.save_pdb_lines(prepared_lines, (PREPARED_PDB_DIR / pdb_file.stem).with_suffix('.pdb'))

        annotation, annotated_lines = volumizer.annotate_pdb_structure(prepared_structure)
        final_pdb_lines = pdb.merge_pdb_lines([prepared_lines, ["END"], annotated_lines])
        annotation_df = utils.make_annotation_dataframe(annotation)

        utils.save_annotation_dataframe(annotation_df, (ANNOTATED_DF_DIR / pdb_file.stem).with_suffix('.json'))
        print(f"Annotation saved to dataframe as: {(ANNOTATED_DF_DIR / pdb_file.stem).with_suffix('.json')}")
        pdb.save_pdb_lines(final_pdb_lines, (ANNOTATED_PDB_DIR / pdb_file.stem).with_suffix('.pdb'))
        print(f"Annotation saved to PDB as: {(ANNOTATED_PDB_DIR / pdb_file.stem).with_suffix('.pdb')}")
        print(f"Quick annotation output:")
        print(annotation_df)


def main(
    volumize_input: str = typer.Argument(
        ..., help="PDB ID, PDB file, file with one PDB ID per line, or folder containing PDB files"
    ),
    resolution: float = typer.Option(VOXEL_SIZE, help="Edge-length of voxels used to discretize the structure."),
    non_protein: bool = typer.Option(False, help="Include non-protein residues during the calculations."),
    jobs: int = typer.Option(1, help="Number of threads to use."),
):
    """
    Find pores and cavities in the supplied PDB files.
    """
    utils.setup_dirs()
    utils.set_resolution(resolution)
    utils.set_non_protein(non_protein)
    utils.add_protein_components(rcsb.load_components())
    input_type = cli.guess_input_type(volumize_input)

    if input_type == "pdb_id":
        volumize_pdb_id(volumize_input)
    elif input_type == "pdb_file":
        pdb_file = Path(volumize_input)
        volumize_pdb_file(pdb_file)
    elif input_type == "id_file":
        with open(volumize_input, mode="r", encoding="utf-8") as id_file:
            pdb_ids = [line.strip() for line in id_file.readlines()]
        tasks = [[pdb_id] for pdb_id in pdb_ids]
        with multiprocessing.Pool(processes=jobs) as pool:
            pool.starmap(volumize_pdb_id, tasks)
    elif input_type == "pdb_dir":
        pdb_files = Path(volumize_input).glob("*.pdb")
        tasks = [[pdb_file] for pdb_file in pdb_files]
        with multiprocessing.Pool(processes=jobs) as pool:
            pool.starmap(volumize_pdb_file, tasks)
    else:
        raise RuntimeError("File mode not implemented")


if "__main__" in __name__:
    typer.run(main)
