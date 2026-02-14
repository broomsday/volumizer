#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if test -f "voxel.so"; then
    rm voxel.so
fi
if test -f "fib_sphere.so"; then
    rm fib_sphere.so
fi

cc -fPIC -shared -o voxel.so voxel.c
cc -fPIC -shared -o fib_sphere.so fib_sphere.c
