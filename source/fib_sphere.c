// fib_sphere.c


#include <math.h>
#include <stdio.h>


float GOLDEN_RATIO = (1.0 + sqrt(5.0)) / 4.0;
float PI = 3.14159265358979323846;


float *fibonacci_sphere(float radius, float x, float y, float z, int samples, float* coords) {
    for (int i = 0; i < samples; i++) {
        float phi = acos(1 - 2 * (i + 0.5) / samples);
        float theta = (PI * i) / GOLDEN_RATIO;

        coords[i] = x + (cos(theta) * sin(phi) * radius);
        coords[samples + i] = y + (sin(theta) * sin(phi) * radius);
        coords[2 * samples + i] = z + (cos(phi) * radius); 
    }

    return coords;
}
