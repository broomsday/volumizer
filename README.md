`volumizer` discovers and annotates occluded volumes in proteins including:
- `cavities`: Volumes within a protein that do not make any contacts with bulk solvent. Useful for e.g. carrying cargo.
- `pockets`: Volumes on the protein surface that make a single contact with bulk solvent.  Useful for e.g. ligand binding or catalysis.
- `pores`: Volumes connecting two bulk solvent surfaces.  Useful for e.g. filtering solutes.
- `hubs`: Volumes connecting more than two bulk solvent surfaces.  Useful for e.g. containing a small reaction volume.

# Example Identified Volume

Here is shown an example pore identified and annotated (red) in PDB 4JPN (green).
![image](images/pore_annotation.png)

The same pore volume annotated (red) shown as a slice through the protein (grey).
![image](images/pore_slice.png)

The dataframe output shows volume/dimensinos of each occluded volume in the structure
|    | id |  type  |  volume  |    x    |    y    |    z    |
|----|----|--------|----------|---------|---------|---------|
| 0  | 0  |   pore | 38286.0  | 108.221 |  38.574 |  36.310 |
| 1  | 0  | pocket |   189.0  |   8.214 |   5.628 |   0.000 |
| 2  | 1  | pocket |   162.0  |   6.635 |   3.843 |   2.840 |
| 3  | 2  | pocket |   162.0  |   7.298 |   4.002 |   2.557 |
| 4  | 3  | pocket |   162.0  |  10.757 |   2.701 |   0.000 |
| 5  | 4  | pocket |   135.0  |   6.000 |   6.000 |   0.000 |

# Installation
The package is published on PyPi, `pip install volumizer`.

## Developing the Volumizer Package
`biotite==0.37.0` does not support Python 3.14. Use Python 3.10 or 3.11 with `uv`.

If you want to develop the package:
1. `uv python install 3.11`
2. `uv sync --python 3.11 --group test`
3. `bash src/compile_c_libs.sh`
4. `uv run --python 3.11 pytest`

Optional native scaffold (Phase 1):
1. `uv sync --python 3.11 --group test --group native`
2. Install Rust toolchain (`cargo`, `rustc`).
3. `uv run --python 3.11 maturin develop --manifest-path native/Cargo.toml`

Backend selection:
- `VOLUMIZER_BACKEND=python` (default): use Python/C-ctypes implementation
- `VOLUMIZER_BACKEND=auto`: use native if importable, otherwise Python
- `VOLUMIZER_BACKEND=native`: require native module and fail if unavailable

# Usage
Using the test file `tests/pdbs/4jpn.pdb` try out the following:

## Volumize a PDB and Save the Volumized PDB and DataFrame
Performing end-to-end loading, cleaning, volumizing, and saving is done with a single convenience function:

```
from volumizer import volumizer

volumizer.volumize_pdb_and_save("my_input.pdb", "volumized_pdb.pdb", "volumized_df.json")
```

## Load a PDB as a Biotite Structure, Clean, Volumize, and Save Output
If you want access to the individual end-products: volume dataframe, the input structure after cleaning, and the structure of the volumes:

```
from volumizer import volumizer, pdb

pdb_structure = pdb.load_structure("my_input.pdb")
volumes_df, cleaned_structure, volumes_structure = volumizer.volumize_structure(pdb_structure)

# take the cleaned input and annotated volumes and convert them to a PDB format string and then save
#  modify the `deliminator` to suit your visualization preference
#  e.g. the default "END" allows Pymol to load the resulting PDB file as two separate objects, one for the cleaned input, and one for the volumes
pdb_lines = pdb.make_volumized_pdb_lines([cleaned_structure, volumes_structure], deliminator="END")
pdb.save_pdb_lines(pdb_lines, "volumized_pdb.pdb")

volumes_df.to_json("volumized_df.json")
```

## Changing Resolution, and Beyond
If you are interested in additional control over the volumizing method:
- the resolution of the voxels can be changed
- cleaning can be skipped

Note: the default voxel resolution is 3.0 Angstroms, which gives sensible results in the majority of cases.
Higher resolutions especially < 2.0 Angstroms will often find small paths through a protein structure, making e.g. cavities look like pores, etc.
Lower resolutions are faster to compute, but may begin to under-estimate the true volume of solvent occluded elements.

Note: by default all residues that make L- or D- peptide bonds are retained through cleaning (e.g. Non-canonicals are kept, even if they are heteroatoms in PDB structure).
By constrast all non-covalently attached residues are removed.  Currently glycan residues are also removed as they make non-peptide bonds, below is shown an example of how would would retain glycans.

```
from volumizer import volumizer, pdb, utils

utils.set_resolution(2.0)

pdb_structure = pdb.load_structure("my_input.pdb")
cleaned_structure = volumizer.prepare_pdb_structure(pdb_structure)  # skip this if you want to keep the exact input structure
volumes_df, volumes_structure = volumizer.annotate_structure_volumes(cleaned_structure)

# take the cleaned input and annotated volumes and convert them to a PDB format string and then save
#  modify the `deliminator` to suit your visualization preference
#  e.g. the default "END" allows Pymol to load the resulting PDB file as two separate objects, one for the cleaned input, and one for the volumes
pdb_lines = pdb.make_volumized_pdb_lines([cleaned_structure, volumes_structure], deliminator="END")
pdb.save_pdb_lines(pdb_lines, "volumized_pdb.pdb")

volumes_df.to_json("volumized_df.json")
```

# How It Works
`volumizer` identifies hydrated volumes in a protein structure that are not fully solvent exposed, e.g. a binding pocket.
It then computes the volume and dimensions of these and outputs that information along with an annotated 
version of the input PDB showing where these volumes are (which can be visualized in e.g. Pymol, Chimerax, etc.).

## Identifying hydrated volumes
1. A large voxel-grid is built around the atoms of the protein or other structure supplied.
2. All voxels within a van der Waals radius of a protein or other atom is flagged as being `non-solvent`
3. Remaining `solvent` atoms are then broken into two groups: `bulk solvent` and `occluded volumes`
   This is done by tracing a vector along each ordinal axis from a given voxel and if 2 or more of these vectors would
   cross a `non-solvent` voxel, then the query voxel is identified as an `occluded volume` to be further analyzed
   otherwise it is considered `bulk solvent`
   The intention is to identify points on the grid that are outside the protein as `bulk solvent`
4. All `occluded volume` voxels are then grouped into a number of continuous volumes
5. For each continous volume the number of distinct surfaces that contact `bulk solvent` voxels is computed and used
   to indicate the volume type:
   0 surfaces interacting with solvent: `cavity`
   1 surface interacting with solvent: `pocket`
   2 surfaces interacting with solvent: `pore`
   3+ surfaces interacting with solvent: `hub`

## Annotations

### Pandas DataFrame
Annotations are given as a pandas data frame saved as a .json file.  The annotation lists all hydrated volumes ordered by
total volume, giving the type of volume, and dimensions.

### PDB File
The input PDB file will be annotated by adding `atoms` to represent the hydrated volumes.  The ATOM entries
contain several points of information about the volume from which they come:

Type of volume: the residue name encodes the type of volume in 3-letter code
- `OCC` for `occluded`
- `CAV` for `cavity`
- `POK` for `pocket`
- `POR` for `pore`
- `HUB` for `hub`

Surface of the hydrated volume that interacts with bulk solvent: this is indicated by a B-factor of 50.0, whereas
the remainder of the volume (that does not interact with the bulk solvent) has a value of 0.0.

All `atoms` of a particular volume are grouped under the same residue number.
