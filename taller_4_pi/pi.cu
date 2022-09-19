#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <chrono>

double PI = 3.14159265358979323846;

/**
 * Función para obtener los parametros pasados desde consola
 */
double get_flag_value(int argc, char **argv, char flag[])
{
    double n = -1, arg = 0;
    for (int i = 0; i < argc; i++)
    {
        if (strcmp(flag, argv[i]) == 0 && !arg)
        {
            n = atoi(argv[i + 1]);
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

/**
 * Aproxima el valor de π utilizando la formula de Leibniz
 */
__global__ void pi_gpu(double *cuda_values, int N) {
    unsigned i = threadIdx.x + (blockDim.x * blockIdx.x);

    double alternate = i % 2 == 0 ? 1 : -1;

    if (i < N)
    {
        cuda_values[i] = alternate / (2 * i + 1);
    }
}

int main(int argc, char **argv)
{
    char n_flag[] = "-n";

    int n, BLOCK_SIZE = 1024;

    n = get_flag_value(argc, argv, n_flag);

    double *cuda_values, *values;

    cudaMalloc((void **)&cuda_values, n * sizeof(double));
    values = (double *)malloc(n * sizeof(double));

    /**
     * inicializa los eventos
     */
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    double approx = 0, error;

    int grid_size = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

    /**
     * Empieza a medir el tiempo que le toma al programa
     */
    cudaEventRecord(start);

    pi_gpu<<<grid_size, BLOCK_SIZE>>>(cuda_values, n);

    /**
     * Copia los resultados calculados en la GPU (cuda_values) en la CPU (values)
     */
    cudaMemcpy(values, cuda_values, n * sizeof(double), cudaMemcpyDeviceToHost);

    /**
     * Sincroniza los kernels para que todos hayan acabado antes de continuar
     */
    cudaDeviceSynchronize();

    /**
     * Realiza la sumatoria para obtener el valor de pi
     */
    for (int i = 0; i < n; i++) {
        approx += values[i];
    }

    approx *= 4;
    /**
     * Obtiene el error de la aproximación
     */
    error = fabs((approx - PI) / PI);

    /**
     * Termina de medir el tiempo
     */
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    printf("Approximated value = %f (Tiempo empleado = %fs, %d iteraciones)\n", approx, milliseconds / 1000, n);
    printf("Real value = %f\n", PI);
    printf("Error = %f\n", error);

    /**
     * Libera la memoria
     */
    cudaFree(cuda_values);
    free(values);
    return 0;
}