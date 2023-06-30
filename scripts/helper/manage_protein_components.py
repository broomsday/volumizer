"""
Given a text file of components, build a python constants file for these.
"""


import typer

from volumizer import rcsb
from volumizer.paths import PROTEIN_COMPONENTS_FILE, PROTEIN_COMPONENTS_SCRIPT


app = typer.Typer()


@app.command()
def make_components_script():
    """
    Write the RCSB components file into a python script to hold as a set.
    """
    with open(PROTEIN_COMPONENTS_FILE, mode="r", encoding="utf-8") as fi:
        components = set([component.rstrip("\n") for component in fi.readlines()])

    component_out = "PROTEIN_COMPONENTS = {\n"
    component_out += ",\n".join([f'\t"{component}"' for component in components])
    component_out += ",\n}\n"

    with open(PROTEIN_COMPONENTS_SCRIPT, mode="w", encoding="utf-8") as fo:
        fo.write(component_out)


@app.command()
def update_components_file():
    """
    Freshly download all the components and update the components file.
    """
    rcsb.get_components()
    make_components_script()


if __name__ == "__main__":
    app()