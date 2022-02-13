"""
Functions to manipulate and analyze voxels.
"""


from copy import deepcopy
import ctypes

import numpy as np
import pandas as pd

from tqdm import tqdm
from pyntcloud import PyntCloud
from pyntcloud.structures.voxelgrid import VoxelGrid

from pore import utils, align
from pore.constants import OCCLUDED_DIMENSION_LIMIT, MIN_NUM_VOXELS, DIAGONAL_NEIGHBORS
from pore.types import VoxelGroup
from pore.paths import C_CODE_DIR


if utils.using_performant():
    VOXEL_C_PATH = C_CODE_DIR / "voxel.so"
    VOXEL_C = ctypes.CDLL(str(VOXEL_C_PATH.absolute()))


def coords_to_point_cloud(coords: pd.DataFrame) -> PyntCloud:
    """
    Produce a point-cloud from atomic coordinates.
    """
    return PyntCloud(coords)


def add_voxel_grid(cloud: PyntCloud) -> tuple[PyntCloud, str]:
    """
    Generate a voxel grid surrounding the point-cloud.
    """
    voxel_grid_id = cloud.add_structure(
        "voxelgrid",
        size_x=utils.VOXEL_SIZE,
        size_y=utils.VOXEL_SIZE,
        size_z=utils.VOXEL_SIZE,
        regular_bounding_box=False,
    )
    return cloud, voxel_grid_id


def get_voxel_grid(cloud: PyntCloud, voxel_grid_id: str) -> VoxelGrid:
    """
    Generate an array representing a binary voxel grid where zero represents no protein
    atoms in that voxel (e.g. solvent) and a non-zero represents a protein atom in that voxel.
    """
    return cloud.structures[voxel_grid_id]


def get_protein_solvent_voxel_array(voxel_grid: VoxelGrid) -> np.ndarray:
    """
    Generate a 3D array of x[y[z]] as the voxel coordinate and the value of that array entry
    being zero for no protein and one when a protein atom is in that voxel.
    """
    return voxel_grid.get_feature_vector(mode="binary")


def get_protein_and_solvent_voxels(
    binary_voxel_array: np.ndarray,
    voxel_grid_dimensions: np.ndarray,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    From the voxel grid, return the indices of the protein containing voxels,
    and those not containing protein (e.g. solvent).
    """
    protein_voxels = np.nonzero(binary_voxel_array)
    solvent_voxels = np.nonzero(binary_voxel_array == 0)

    protein_voxel_indices = compute_voxel_indices(protein_voxels, voxel_grid_dimensions)
    solvent_voxel_indices = compute_voxel_indices(solvent_voxels, voxel_grid_dimensions)

    return (
        VoxelGroup(
            voxels=protein_voxels,
            indices=protein_voxel_indices,
            num_voxels=len(protein_voxel_indices),
            voxel_type="protein",
            volume=compute_voxel_group_volume(len(protein_voxel_indices)),
        ),
        VoxelGroup(
            voxels=solvent_voxels,
            indices=solvent_voxel_indices,
            num_voxels=len(solvent_voxel_indices),
            voxel_type="solvent",
            volume=compute_voxel_group_volume(len(solvent_voxel_indices)),
        ),
    )


def get_occluded_dimensions(
    query_voxel: tuple[int, int, int],
    possible_occluding_x: list[int],
    possible_occluding_y: list[int],
    possible_occluding_z: list[int],
) -> list[int]:
    """
    Determine how many ordinal axes are occluded.
    """
    occluded_dimensions = [0, 0, 0, 0, 0, 0]

    if len(possible_occluding_z) != 0:
        if min(possible_occluding_z) < query_voxel[2]:
            occluded_dimensions[4] = 1
        if max(possible_occluding_z) > query_voxel[2]:
            occluded_dimensions[5] = 1
    if len(possible_occluding_y) != 0:
        if min(possible_occluding_y) < query_voxel[1]:
            occluded_dimensions[2] = 1
        if max(possible_occluding_y) > query_voxel[1]:
            occluded_dimensions[3] = 1
    if len(possible_occluding_x) != 0:
        if min(possible_occluding_x) < query_voxel[0]:
            occluded_dimensions[0] = 1
        if max(possible_occluding_x) > query_voxel[0]:
            occluded_dimensions[1] = 1

    return occluded_dimensions


def is_buried(occluded_dimensions: list[int]) -> bool:
    """
    If 5 or 6 dimensions are occluded, return True.
    If less than 4 dimensions are occluded, return False.
    If exactly 4 dimensions are occluded, return True if the two unoccluded dimensions are the same axis,
    False otherwise
    """

    if occluded_dimensions.count(1) > OCCLUDED_DIMENSION_LIMIT:
        return True
    elif occluded_dimensions.count(1) == OCCLUDED_DIMENSION_LIMIT:
        if occluded_dimensions[0] == 0 and occluded_dimensions[1] == 0:
            return True
        elif occluded_dimensions[2] == 0 and occluded_dimensions[3] == 0:
            return True
        elif occluded_dimensions[4] == 0 and occluded_dimensions[5] == 0:
            return True

    return False


def build_planar_voxel_coordinate_arrays(
    voxels: tuple[np.ndarray, ...], voxel_grid_dimensions: np.ndarray
) -> tuple[list[list[list[int]]], list[list[list[int]]], list[list[list[int]]]]:
    """
    For each two-dimension pair, construct a list of coordinates in the 3rd dimension that match the first two dimensions.
    """
    # TODO collapse down the final lists to min() and max() only as two different entries and use that for comparison in the downstream function
    z_array = [[list() for y in range(voxel_grid_dimensions[1])] for x in range(voxel_grid_dimensions[0])]
    y_array = [[list() for z in range(voxel_grid_dimensions[2])] for x in range(voxel_grid_dimensions[0])]
    x_array = [[list() for z in range(voxel_grid_dimensions[2])] for y in range(voxel_grid_dimensions[1])]

    for i in range(voxels[0].size):
        z_array[voxels[0][i]][voxels[1][i]].append(voxels[2][i])
        y_array[voxels[0][i]][voxels[2][i]].append(voxels[1][i])
        x_array[voxels[1][i]][voxels[2][i]].append(voxels[0][i])

    return z_array, y_array, x_array


def get_exposed_and_buried_voxels(
    solvent_voxels: VoxelGroup,
    protein_voxels: VoxelGroup,
    voxel_grid_dimensions: np.ndarray,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    Use simple geometric heuristics to determine if a group of solvent voxels is buried or exposed.
    """
    buried_voxels = ([], [], [])
    exposed_voxels = ([], [], [])

    z_array, y_array, x_array = build_planar_voxel_coordinate_arrays(protein_voxels.voxels, voxel_grid_dimensions)

    for i in tqdm(range(solvent_voxels.voxels[0].size), desc="Searching voxels"):
        query_voxel = (solvent_voxels.voxels[0][i], solvent_voxels.voxels[1][i], solvent_voxels.voxels[2][i])

        occluded_dimensions = get_occluded_dimensions(
            query_voxel,
            x_array[query_voxel[1]][query_voxel[2]],
            y_array[query_voxel[0]][query_voxel[2]],
            z_array[query_voxel[0]][query_voxel[1]],
        )

        if is_buried(occluded_dimensions):
            buried_voxels[0].append(query_voxel[0])
            buried_voxels[1].append(query_voxel[1])
            buried_voxels[2].append(query_voxel[2])
        else:
            exposed_voxels[0].append(query_voxel[0])
            exposed_voxels[1].append(query_voxel[1])
            exposed_voxels[2].append(query_voxel[2])

    buried_voxel_indices = compute_voxel_indices(buried_voxels, voxel_grid_dimensions)
    exposed_voxel_indices = compute_voxel_indices(exposed_voxels, voxel_grid_dimensions)

    return (
        VoxelGroup(
            voxels=(np.array(exposed_voxels[0]), np.array(exposed_voxels[1]), np.array(exposed_voxels[2])),
            indices=exposed_voxel_indices,
            num_voxels=len(exposed_voxel_indices),
            voxel_type="exposed",
            volume=compute_voxel_group_volume(len(exposed_voxel_indices)),
        ),
        VoxelGroup(
            voxels=(np.array(buried_voxels[0]), np.array(buried_voxels[1]), np.array(buried_voxels[2])),
            indices=buried_voxel_indices,
            num_voxels=len(buried_voxel_indices),
            voxel_type="buried",
            volume=compute_voxel_group_volume(len(buried_voxel_indices)),
        ),
    )


def get_first_shell_exposed_voxels(
    exposed_voxels: VoxelGroup, buried_voxels: VoxelGroup, voxel_grid: VoxelGrid
) -> VoxelGroup:
    """
    Agglomerate buried solvent voxels into pores, pockets, cavities, and simply occluded.
    """
    # TODO

    return VoxelGroup()


def is_neighbor_voxel(
    voxel_one: tuple[np.int64, ...], voxel_two: tuple[np.int64, ...], diagonal_neighbors: bool = DIAGONAL_NEIGHBORS
) -> bool:
    """
    Given two voxels return True if they are ordinal neighbors.
    """
    # below approach is ~2x faster than using a list comprehension, presumably because of early exit
    sum_differences = 0
    max_difference = 0
    for dimension in range(3):
        difference = abs(voxel_one[dimension] - voxel_two[dimension])

        if difference > 1:
            return False

        sum_differences += difference
        if difference > max_difference:
            max_difference = difference

    # a voxel is an ordinal neighbor when the sum of the absolute differences in axes indices is 1
    if sum_differences == 1:
        return True
    # alternatively a voxel is a diagonal neighbour when the max distance on any axis is 1
    # here we only count being diagonal on a plane, not in all 3 dimensions
    #   which would be a longer distance between voxels
    elif (diagonal_neighbors) and (sum_differences == 2) and (max_difference == 1):
        return True

    return False


def get_single_voxel(voxels: tuple[np.ndarray, ...], index: int) -> tuple[np.int64, ...]:
    """
    Given a set of voxels return just one voxel at index.
    """
    return (
        voxels[0][index],
        voxels[1][index],
        voxels[2][index],
    )


def compute_voxel_indices(voxels: tuple[np.ndarray, ...], grid_dimensions: np.ndarray) -> set[int]:
    """
    Given a 3D array of voxels, compute the 1D set of their indices.
    """
    # TODO: is this built-in easier and correct?: x, y, z = np.unravel_index(voxel, self.x_y_z)
    return {
        (voxels[0][i] * grid_dimensions[1] * grid_dimensions[2] + voxels[1][i] * grid_dimensions[2] + voxels[2][i])
        for i in range(len(voxels[0]))
    }


def compute_voxel_group_volume(num_voxels: int) -> float:
    """
    Return the volume of a number of voxels
    """
    return num_voxels * utils.VOXEL_VOLUME


def breadth_first_search_python(voxels: tuple[np.ndarray, ...], searchable_indices: set[int]) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    searchable_indices = deepcopy(searchable_indices)
    start_index = searchable_indices.pop()
    queue_indices = set([start_index])
    neighbor_indices = set([start_index])

    while len(queue_indices) > 0:
        current_index = queue_indices.pop()
        current_voxel = get_single_voxel(voxels, current_index)
        for searched_index in searchable_indices:
            searched_voxel = get_single_voxel(voxels, searched_index)
            if is_neighbor_voxel(current_voxel, searched_voxel):
                queue_indices.add(searched_index)
                neighbor_indices.add(searched_index)
        searchable_indices -= neighbor_indices

    return neighbor_indices


def breadth_first_search_c(voxels: tuple[np.ndarray, ...], searchable_indices: set[int]) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    c_searchable_indices = list(deepcopy(searchable_indices))

    # make the voxels compatible with C
    num_voxels = len(voxels[0])
    voxels_x = (ctypes.c_int * num_voxels)(*voxels[0])
    voxels_y = (ctypes.c_int * num_voxels)(*voxels[1])
    voxels_z = (ctypes.c_int * num_voxels)(*voxels[2])

    # generate other inputs needed for C function
    c_searchable_indices = (ctypes.c_int * len(c_searchable_indices))(*c_searchable_indices)
    initial_indices = (ctypes.c_int * num_voxels)(*([-1] * num_voxels))

    # run the breadth-first search
    VOXEL_C.breadth_first_search.restype = ctypes.POINTER(ctypes.c_int * num_voxels)
    neighbor_indices_c = VOXEL_C.breadth_first_search(
        ctypes.byref(voxels_x),
        ctypes.byref(voxels_y),
        ctypes.byref(voxels_z),
        ctypes.byref(c_searchable_indices),
        len(searchable_indices),
        ctypes.byref(initial_indices),
    )

    return set([i for i in neighbor_indices_c.contents if i != -1])


def breadth_first_search(
    voxels: tuple[np.ndarray, ...], searchable_indices: set[int], performant: bool = utils.using_performant()
) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    if performant:
        return breadth_first_search_c(voxels, searchable_indices)

    return breadth_first_search_python(voxels, searchable_indices)


def get_agglomerated_type(
    query_indices: set[int], buried_voxels: tuple[np.ndarray, ...], exposed_voxels: tuple[np.ndarray, ...]
) -> tuple[set[int], str]:
    """
    Find "surface" voxels, being buried voxel in direct contact with an exposed voxel. Three possibilites:

    a) there are no "surface" voxels -> this is a cavity
    b) all "surface" voxels can be agglomerated (BFS) into a single group -> this is a pocket
    c) the "surafce" voxels cannot be agglomerated (BFS) into just one group -> this is a pore
    """
    direct_surface_indices = set()
    for query_index in query_indices:
        query_voxel = get_single_voxel(buried_voxels, query_index)
        for exposed_voxel_index in range(len(exposed_voxels[0])):
            exposed_voxel = get_single_voxel(exposed_voxels, exposed_voxel_index)
            if is_neighbor_voxel(query_voxel, exposed_voxel):
                direct_surface_indices.add(query_index)
                break

    # if there are no surface contacts, this must be a cavity
    if len(direct_surface_indices) == 0:
        return direct_surface_indices, "cavity"

    # add all voxels that are neighbours to the direct surface voxels
    neighbor_surface_indices = set()
    for query_index in query_indices - direct_surface_indices:
        query_voxel = get_single_voxel(buried_voxels, query_index)
        for surface_index in direct_surface_indices:
            surface_voxel = get_single_voxel(buried_voxels, surface_index)
            if is_neighbor_voxel(surface_voxel, query_voxel):
                neighbor_surface_indices.add(query_index)
                break

    # we pass all surface indices and the buried voxels for BFS
    # If there is more than one distinct surface (e.g. a pore) then
    #   BFS will terminate after agglomerating just one of the surfaces
    #   and `single_surface_indices` will be less than `surface_indices`
    surface_indices = direct_surface_indices.union(neighbor_surface_indices)
    single_surface_indices = breadth_first_search(buried_voxels, surface_indices)
    if len(single_surface_indices) < len(surface_indices):
        return direct_surface_indices, "pore"

    return direct_surface_indices, "pocket"


def get_voxel_group_center(voxel_indices: set[int], voxel_grid: VoxelGrid) -> np.ndarray:
    """
    Compute the center of geometry for the voxel group.
    """
    voxel_coords = np.array([voxel_grid.voxel_centers[voxel_index] for voxel_index in voxel_indices])
    return np.mean(voxel_coords, axis=0)


def get_voxel_group_axial_lengths(voxel_indices: set[int], voxel_grid: VoxelGrid) -> list[float]:
    """
    Align the voxel group to it's principal axes, then compute the maximum length along each axis.
    """
    voxel_coords = np.array([voxel_grid.voxel_centers[voxel_index] for voxel_index in voxel_indices])

    rotation, translation = align.get_principal_axis_alignment_translation_rotation(voxel_coords)

    # apply the rotation and translation
    voxel_coords = align.rotate_and_translate_coords(voxel_coords, rotation, translation)

    # compute the maximum length across each axis
    xs = [voxel_coord[0] for voxel_coord in voxel_coords]
    ys = [voxel_coord[1] for voxel_coord in voxel_coords]
    zs = [voxel_coord[2] for voxel_coord in voxel_coords]

    return sorted([max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)], reverse=True)


def get_pores_pockets_cavities_occluded(
    buried_voxels: VoxelGroup, exposed_voxels: VoxelGroup, voxel_grid: VoxelGrid
) -> tuple[dict[int, VoxelGroup], dict[int, VoxelGroup], dict[int, VoxelGroup], dict[int, VoxelGroup]]:
    """
    Agglomerate buried solvent voxels into pores, pockets, cavities, and simply occluded.
    """
    buried_indices = set(range(buried_voxels.voxels[0].size))
    agglomerated_indices = set()

    pores, pockets, cavities, occluded = {}, {}, {}, {}
    pore_id, pocket_id, cavity_id, occluded_id = 0, 0, 0, 0

    while len(agglomerated_indices) < len(buried_indices):
        remaining_indices = buried_indices - agglomerated_indices

        # perform BFS over the remaining indices
        agglomerable_indices = breadth_first_search(buried_voxels.voxels, remaining_indices)
        # iterate our counter of finished indices
        agglomerated_indices = agglomerated_indices.union(agglomerable_indices)

        # if too small, don't assign direct surface indices and assign type "occluded"
        if len(agglomerable_indices) <= MIN_NUM_VOXELS:
            direct_surface_indices = set()
            agglomerated_type = "occluded"
        else:
            # identify what these agglomerated voxels are
            direct_surface_indices, agglomerated_type = get_agglomerated_type(
                agglomerable_indices, buried_voxels.voxels, exposed_voxels.voxels
            )

        # get the surface voxels for use in getting their voxel-grid indices
        surface_voxels = (
            np.array([buried_voxels.voxels[0][index] for index in direct_surface_indices]),
            np.array([buried_voxels.voxels[1][index] for index in direct_surface_indices]),
            np.array([buried_voxels.voxels[2][index] for index in direct_surface_indices]),
        )

        # get the voxels
        volume_voxels = (
            np.array([buried_voxels.voxels[0][index] for index in agglomerable_indices]),
            np.array([buried_voxels.voxels[1][index] for index in agglomerable_indices]),
            np.array([buried_voxels.voxels[2][index] for index in agglomerable_indices]),
        )
        # get the voxel indices
        volume_indices = compute_voxel_indices(volume_voxels, voxel_grid.x_y_z)
        # create the voxelgroup
        voxel_group = VoxelGroup(
            voxels=volume_voxels,
            indices=volume_indices,
            surface_indices=compute_voxel_indices(surface_voxels, voxel_grid.x_y_z),
            num_voxels=len(volume_indices),
            voxel_type=agglomerated_type,
            volume=compute_voxel_group_volume(len(volume_indices)),
            center=get_voxel_group_center(volume_indices, voxel_grid),
            axial_lengths=get_voxel_group_axial_lengths(volume_indices, voxel_grid),
        )

        # add the voxelgroup depending on type and increment the type counter
        if agglomerated_type == "cavity":
            cavities[cavity_id] = voxel_group
            cavity_id += 1
        elif agglomerated_type == "pore":
            pores[pore_id] = voxel_group
            pore_id += 1
        elif agglomerated_type == "pocket":
            pockets[pocket_id] = voxel_group
            pocket_id += 1
        elif agglomerated_type == "occluded":
            occluded[occluded_id] = voxel_group
            occluded_id += 1

    return pores, pockets, cavities, occluded
