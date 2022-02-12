rm source/voxel_compute.so
rm source/fib_sphere.so
cc -fPIC -shared -o source/voxel_compute.so source/voxel_compute.c
cc -fPIC -shared -o source/fib_sphere.so source/fib_sphere.c