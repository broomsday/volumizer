// voxel_compute.c


#include <stdlib.h>


int *get_single_voxel(int* voxels_x, int* voxels_y, int* voxels_z, int index, int* voxel_indices) {
    voxel_indices[0] = voxels_x[index];
    voxel_indices[1] = voxels_y[index];
    voxel_indices[2] = voxels_z[index];

    return voxel_indices;
}


int is_neighbor_voxel(int* voxel_one, int* voxel_two, int diagonal_neighbors) {
    int sum_differences = 0;
    int max_difference = 0;
    for (int i = 0; i < 3; i++) {
        int difference = abs(voxel_one[i] - voxel_two[i]);
        sum_differences += difference;
        if (difference > max_difference) {
            max_difference = difference;
        }
    }

    if (sum_differences == 1) {
        return 1;
    }
    else if ((diagonal_neighbors == 1) && (sum_differences == 2) && (max_difference == 1)) {
        return 1;
    }
    return 0;
}


int count_positive_ints(int values[]) {
    int num_positive = 0;
    for (int i = 0; i < (sizeof(*values) / sizeof(int)); i++) {
        if (values[i] >= 0) {
            num_positive += 1;
        }
    }

    return num_positive;
}


int *breadth_first_search(int* voxels_x, int* voxels_y, int* voxels_z, int* input_searchable_indices, int num_searchable_indices, int* neighbor_indices) {
    int current_voxel_indices[3];
    int searched_voxel_indices[3];

    int searchable_indices[num_searchable_indices];
    int queue_indices[num_searchable_indices];

    for (int i = 0; i < num_searchable_indices; i++) {
        searchable_indices[i] = input_searchable_indices[i];
        queue_indices[i] = -1;
    }

    int queue_position = 0;
    queue_indices[queue_position] = searchable_indices[0];
    neighbor_indices[0] = searchable_indices[0];
    searchable_indices[0] = -1;

    while (count_positive_ints(queue_indices) > 0) {
        int current_index = queue_indices[queue_position];
        queue_indices[queue_position] = -1;
        queue_position -= 1;

        int* current_voxel = get_single_voxel(voxels_x, voxels_y, voxels_z, current_index, current_voxel_indices);
        for (int i = 0; i < num_searchable_indices; i++) {
            if (searchable_indices[i] == -1) {
                continue;
            }
            int* searched_voxel = get_single_voxel(voxels_x, voxels_y, voxels_z, searchable_indices[i], searched_voxel_indices);
            if (is_neighbor_voxel(current_voxel, searched_voxel, 1)) {
                queue_position += 1;
                queue_indices[queue_position] = searchable_indices[i];
                neighbor_indices[i] = searchable_indices[i];
            }
        }
        for (int i = 0; i < num_searchable_indices; i++) {
            if (neighbor_indices[i] != -1) {
                searchable_indices[i] = -1;
            }
        }
    }

    return neighbor_indices;
}
