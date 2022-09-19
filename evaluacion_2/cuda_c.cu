% % cuda-- name test_hw.cu
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <string.h>
#include <cuda.h>
#include <curand.h>

#define BLOCK_SIZE 1024;
#define NUMBER_THREADS_PER_BLOCK 16;

    double
    get_flag_value(int argc, char **argv, char flag[])
{
    double n = -1, arg = 0;
    for (int i = 0; i < argc; i++)
    {
        if (strcmp(flag, argv[i]) == 0 && !arg)
        {
            n = atof(argv[i + 1]);
            arg = 1;
        }
        else
        {
            continue;
        }
    }

    if (n == -1)
    {
        printf("Missing flag.\n");
    }

    return n;
}

double f(double x)
{
    return pow(x, 2);
}

double F(double x)
{
    return pow(x, 3) / 3;
}

__global__ void montecarlo_gpu(double *estimations, double a, double b, int N)
{
    unsigned idx = threadIdx.x + (blockDim.x * blockIdx.x);
    unsigned idy = threadIdx.y + (blockDim.y * blockIdx.y);
    int i = idx * N + idy;
    if (i < N)
    {
        estimations[i] = (b - a) * f(a + (b - a) * 0.2);
    }
}

int main(int argc, char **argv)
{
    char a_flag[] = "-a";
    char b_flag[] = "-b";
    char n_flag[] = "-n";

    double a, b;
    int n;

    a = get_flag_value(argc, argv, a_flag);
    b = get_flag_value(argc, argv, b_flag);
    n = floor(get_flag_value(argc, argv, n_flag));

    double random_estimations[n];

    double *cuda_random_estimations;

    cudaMalloc(&cuda_random_estimations, n * sizeof(int));

    printf("a = %f, b = %f, n = %f\n", a, b, n);
    double approx, real, error;
    dim3 block(NUMBER_THREADS_PER_BLOCK, NUMBER_THREADS_PER_BLOCK);
    dim3 grid(((n + NUMBER_THREADS_PER_BLOCK - 1) / NUMBER_THREADS_PER_BLOCK), ((n + NUMBER_THREADS_PER_BLOCK - 1) / NUMBER_THREADS_PER_BLOCK));

    /**
     * Mide el tiempo de ejecución de la aproximación
     */
    montecarlo_gpu<<<block, threads_per_block>>>(cuda_random_estimations, a, b, n);

    /**
     * calcula la integral de la manera estándar
     */
    real = F(b) - F(a);
    /**
     * Obtiene el error de la aproximación
     */
    error = fabs((approx - real) / real);

    printf("Approximated value = %f (Tiempo empleado = %fs, %f iteraciones)\n", approx, time_spent, n);
    printf("Real value = %f\n", real);
    printf("Error = %f\n", error);

    return 0;
}

dim3 block(NUMBER_THREADS_PER_BLOCK, NUMBER_THREADS_PER_BLOCK);
dim3 grid(((n + NUMBER_THREADS_PER_BLOCK - 1) / NUMBER_THREADS_PER_BLOCK), ((n + NUMBER_THREADS_PER_BLOCK - 1) / NUMBER_THREADS_PER_BLOCK));