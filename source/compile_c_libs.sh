if test -f "source/voxel.so"; then
    rm source/voxel.so
fi
if test -f "source/fib_sphere.so"; then
    rm source/fib_sphere.so
fi
cc -fPIC -shared -o source/voxel.so source/voxel.c
cc -fPIC -shared -o source/fib_sphere.so source/fib_sphere.c