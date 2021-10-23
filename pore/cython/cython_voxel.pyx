import numpy as np

from tqdm import tqdm

OCCLUDED_DIMENSION_LIMIT = 4

def get_num_occluded_dimensions(tuple query_voxel, list possible_occluding_x, list possible_occluding_y, list possible_occluding_z):
    """
    Determine how many ordinal axes are occluded.
    """
    cdef int num_occluded_dimensions = 0

    if len(possible_occluding_z) != 0:
        if min(possible_occluding_z) < query_voxel[0]:
            num_occluded_dimensions += 1
        elif max(possible_occluding_z) > query_voxel[0]:
            num_occluded_dimensions += 1
    if len(possible_occluding_y) != 0:
        if min(possible_occluding_y) < query_voxel[0]:
            num_occluded_dimensions += 1
        elif max(possible_occluding_y) > query_voxel[0]:
            num_occluded_dimensions += 1
    if len(possible_occluding_x) != 0:
        if min(possible_occluding_x) < query_voxel[0]:
            num_occluded_dimensions += 1
        elif max(possible_occluding_x) > query_voxel[0]:
            num_occluded_dimensions += 1

    return num_occluded_dimensions


def get_exposed_and_buried_solvent_voxels(tuple solvent_voxels, tuple protein_voxels):
    """
    Use simple geometric heuristics to determine if a given solvent voxel is buried or exposed.
    """

    cdef tuple buried_solvent_voxels = ([], [], [])
    cdef tuple exposed_solvent_voxels = ([], [], [])

    cdef tuple previous_voxel = (-1, -1, -1)
    for i in tqdm(range(solvent_voxels[0].size), desc="Searching voxels"):
        # TODO numpy.where() is somewhat slow, might be a way to improve
        query_voxel = (solvent_voxels[0][i], solvent_voxels[1][i], solvent_voxels[2][i])
        if query_voxel[0] != previous_voxel[0]:
            match_x_protein_voxel_indices = set(np.where(protein_voxels[0] == query_voxel[0])[0])
        if query_voxel[1] != previous_voxel[1]:
            match_y_protein_voxel_indices = set(np.where(protein_voxels[1] == query_voxel[1])[0])
        if query_voxel[2] != previous_voxel[2]:
            match_z_protein_voxel_indices = set(np.where(protein_voxels[2] == query_voxel[2])[0])
        previous_voxel = query_voxel

        possible_occluding_z_indices = match_x_protein_voxel_indices.intersection(match_y_protein_voxel_indices)
        possible_occluding_y_indices = match_x_protein_voxel_indices.intersection(match_z_protein_voxel_indices)
        possible_occluding_x_indices = match_y_protein_voxel_indices.intersection(match_z_protein_voxel_indices)

        possible_occluding_z = [protein_voxels[2][i] for i in possible_occluding_z_indices]
        possible_occluding_y = [protein_voxels[1][i] for i in possible_occluding_y_indices]
        possible_occluding_x = [protein_voxels[0][i] for i in possible_occluding_x_indices]

        num_occluded_dimensions = get_num_occluded_dimensions(
            query_voxel, possible_occluding_x, possible_occluding_y, possible_occluding_z
        )

        if num_occluded_dimensions >= OCCLUDED_DIMENSION_LIMIT:
            buried_solvent_voxels[0].append(query_voxel[0])
            buried_solvent_voxels[1].append(query_voxel[1])
            buried_solvent_voxels[2].append(query_voxel[2])
        else:
            exposed_solvent_voxels[0].append(query_voxel[0])
            exposed_solvent_voxels[1].append(query_voxel[1])
            exposed_solvent_voxels[2].append(query_voxel[2])

    return (
        (
            np.array(exposed_solvent_voxels[0]),
            np.array(exposed_solvent_voxels[1]),
            np.array(exposed_solvent_voxels[2]),
        ),
        (
            np.array(buried_solvent_voxels[0]),
            np.array(buried_solvent_voxels[1]),
            np.array(buried_solvent_voxels[2]),
        ),
    )
