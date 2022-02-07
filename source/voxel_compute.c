// voxel_compute.c


int *get_single_voxel(int* voxels_x, int* voxels_y, int* voxels_z, int index, int* voxel_indices) {
    voxel_indices[0] = voxels_x[index];
    voxel_indices[1] = voxels_y[index];
    voxel_indices[2] = voxels_z[index];

    return voxel_indices;
}


bool is_neighbor_voxel() {
    
}


/*
// python code for is_neighbor_voxel
def is_neighbor_voxel(voxel_one, voxel_two, diagonal_neighbors: bool = DIAGONAL_NEIGHBORS) -> bool:
    """
    Given two voxels return True if they are ordinal neighbors.
    """
    differences = []
    for dimension in range(3):
        differences.append(abs(voxel_one[dimension] - voxel_two[dimension]))

    # a voxel is an ordinal neighbor when the sum of the absolute differences in axes indices is 1
    if sum(differences) == 1:
        return True
    # alternatively a voxel is a diagonal neighbour when the max distance on any axis is 1
    # here we only count being diagonal on a plane, not in all 3 dimensions
    #   which would be a longer distance between voxels
    elif (diagonal_neighbors) and (sum(differences) == 2) and (max(differences) == 1):
        return True

    return False
*/


int *breadth_first_search(int* indices) {
    indices[0] = 0;
    indices[1] = 19;
    indices[2] = 3;
    indices[3] = 14;

    return indices;
}
