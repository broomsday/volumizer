import time

import typer

from volumizer import voxel, pdb, fib_sphere, utils, paths


def main(code: str = "python", pdb_name: str = "pore"):
    """
    Do performance profiling of pore.voxel.breadth_first_search and related functions
    """
    # load some test data from a PDB of interest
    start_load_time = time.time()
    structure = pdb.load_structure(paths.DATA_DIR / "test_data" / f"{pdb_name}.pdb")
    coords = pdb.get_structure_coords(structure)
    end_load_time = time.time()

    start_fib_time = time.time()
    if code == "python":
        coords = fib_sphere.add_extra_points_python(coords, utils.VOXEL_SIZE)
    elif code == "C":
        coords = fib_sphere.add_extra_points_c(coords, utils.VOXEL_SIZE)
    else:
        raise RuntimeError("Language not implemented")
    end_fib_time = time.time()

    start_voxelize_time = time.time()
    cloud = voxel.coords_to_point_cloud(coords)
    cloud, voxel_grid_id = voxel.add_voxel_grid(cloud)
    voxel_grid = voxel.get_voxel_grid(cloud, voxel_grid_id)
    end_voxelize_time = time.time()

    start_solvent_time = time.time()
    protein_solvent_voxels = voxel.get_protein_solvent_voxel_array(voxel_grid)
    protein_voxels, solvent_voxels = voxel.get_protein_and_solvent_voxels(protein_solvent_voxels, voxel_grid.x_y_z)
    _, buried_voxels = voxel.get_exposed_and_buried_voxels(solvent_voxels, protein_voxels, voxel_grid.x_y_z)
    end_solvent_time = time.time()

    start_bfs_time = time.time()
    buried_indices = set(range(buried_voxels.voxels[0].size))
    agglomerated_indices = set()
    # do BFS until we've agglomerated all indices into neighbour groups
    while len(agglomerated_indices) < len(buried_indices):
        remaining_indices = buried_indices - agglomerated_indices

        # perform BFS over the remaining indices
        if code == "python":
            agglomerable_indices = voxel.breadth_first_search_python(buried_voxels.voxels, remaining_indices)
        elif code == "C":
            agglomerable_indices = voxel.breadth_first_search_c(buried_voxels.voxels, remaining_indices)
        else:
            raise RuntimeError("Language not implemented")

        # iterate our counter of finished indices
        agglomerated_indices = agglomerated_indices.union(agglomerable_indices)
    end_bfs_time = time.time()

    load_time = end_load_time - start_load_time
    fib_time = end_fib_time - start_fib_time
    voxelize_time = end_voxelize_time - start_voxelize_time
    solvent_time = end_solvent_time - start_solvent_time
    bfs_time = end_bfs_time - start_bfs_time

    print(agglomerated_indices)
    print(len(agglomerated_indices))
    print(f"\nLoad time: {load_time}")
    print(f"Fibonacci time: {fib_time}")
    print(f"Preprocess time: {voxelize_time}")
    print(f"Init search time: {solvent_time}")
    print(f"BFS time: {bfs_time}")


if "__main__" in __name__:
    typer.run(main)