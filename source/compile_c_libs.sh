rm source/voxel.so
rm source/fib_sphere.so
cc -fPIC -shared -o source/voxel.so source/voxel.c
cc -fPIC -shared -o source/fib_sphere.so source/fib_sphere.c