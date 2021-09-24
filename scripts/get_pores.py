"""
Download biological assemblies from the RCSB and find those that are pores capable of supporting de novo enzyme design.
"""


import typer

from pore import rcsb, pdb


def main():
    """
    Download biological assemblies, clean and process them, then find pores.
    """

    pdbs_to_download = rcsb.build_pdb_list()
    print(pdbs_to_download)


if "__main__" in __name__:
    typer.run(main)