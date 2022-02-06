// voxel_compute.c


#include <stdlib.h>
#include <stdio.h>


typedef struct VoxelIndices {
    int x;
    int y;
    int z;
} VOXEL_INDICES;

/*
typedef struct Voxels {
    int x[];
    int y[];
    int z[];
} VOXELS;
*/

VOXEL_INDICES *get_single_voxel(int* voxels_x, int* voxels_y, int* voxels_z, int index) {
    VOXEL_INDICES *indices;
    VOXEL_INDICES initial = { voxels_x[index], voxels_y[index], voxels_z[index] };
    indices = malloc(sizeof(VOXEL_INDICES));
    *indices = initial;

    return indices;
}


void free_voxel_indices(VOXEL_INDICES *indices) {
    free(indices);
}