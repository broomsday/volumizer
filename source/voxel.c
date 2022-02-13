// voxel_compute.c


#include <stdlib.h>
#include <stdio.h>


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

        if (difference > 1) {
            return 0;
        }

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
    // initialize arrays to hold return values for `get_single_voxel`
    int current_voxel_indices[3];
    int searched_voxel_indices[3];

    // initialize and set values for a copy of searchable indices and a queue
    int searchable_indices[num_searchable_indices];
    int queue_indices[num_searchable_indices];
    for (int i = 0; i < num_searchable_indices; i++) {
        searchable_indices[i] = input_searchable_indices[i];
        queue_indices[i] = -1;
    }

    // initialize a queue and stack
    int queue_position = 0, indices_position = 0;
    queue_indices[queue_position] = searchable_indices[indices_position];
    neighbor_indices[indices_position] = searchable_indices[indices_position];
    indices_position++;

    while (count_positive_ints(queue_indices) > 0) {
        int current_index = queue_indices[queue_position];
        queue_indices[queue_position] = -1;
        queue_position--;

        int* current_voxel = get_single_voxel(voxels_x, voxels_y, voxels_z, current_index, current_voxel_indices);
        for (int i = indices_position; i < num_searchable_indices; i++) {
            int* searched_voxel = get_single_voxel(voxels_x, voxels_y, voxels_z, searchable_indices[i], searched_voxel_indices);
            if (is_neighbor_voxel(current_voxel, searched_voxel, 1)) {
                queue_position++;;
                queue_indices[queue_position] = searchable_indices[i];

                neighbor_indices[indices_position] = searchable_indices[i];
                searchable_indices[i] = searchable_indices[indices_position];
                indices_position++;
            }
        }
    }

    return neighbor_indices;
}


int *get_neighbor_voxels(int* query_voxels_x, int* query_voxels_y, int* query_voxels_z, int* reference_voxels_x, int* reference_voxels_y, int* reference_voxels_z, int num_query, int num_reference, int* neighbor_indices) {
    // initialize arrays to hold return values for `get_single_voxel`
    int query_voxel_indices[3];
    int reference_voxel_indices[3];

    for (int query_index = 0; query_index < num_query; query_index++) {
        int* query_voxel = get_single_voxel(query_voxels_x, query_voxels_y, query_voxels_z, query_index, query_voxel_indices);
        for (int reference_index = 0; reference_index < num_reference; reference_index++) {
            int* reference_voxel = get_single_voxel(reference_voxels_x, reference_voxels_y, reference_voxels_z, reference_index, reference_voxel_indices);
            if (is_neighbor_voxel(query_voxel, reference_voxel, 1)) {
                neighbor_indices[query_index] = query_index;
                break;
            }
        }
    }

    return neighbor_indices;
}
