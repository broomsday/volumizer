// voxel_compute.c


#include <stdlib.h>
#include <stdio.h>


typedef struct VoxelIndices {
    int a;
    int b;
    int c;
} VOXEL_INDICES;


VOXEL_INDICES *get_single_voxel() {
    VOXEL_INDICES *indices;
    VOXEL_INDICES initial = { 11, 2, 65 };
    indices = malloc(sizeof(VOXEL_INDICES));
    *indices = initial;

    return indices;
}


void free_voxel_indices(VOXEL_INDICES *indices) {
    free(indices);
}