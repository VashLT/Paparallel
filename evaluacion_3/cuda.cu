#include <stdio.h>
#include <stdlib.h>
/**
 * Para utilizar la función pow
 */
#include <math.h>
/**
 * Para parsear los argumentos desde la linea de comandos
 */
#include <string.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <chrono>


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

/**
 * función de la cual se quiere obtener la integral
 */
__device__ double f(double x)
{
    return pow(x, 2);
}

/**
 * Antiderivada de la función original
 */
double F(double x)
{
    return pow(x, 3) / 3;
}

/**
 * Kernel para calcular el método de integración montecarlo
 */
__global__ void montecarlo_gpu(double *estimations, curandState *state, double a, double b, int N)
{
    unsigned i = threadIdx.x + (blockDim.x * blockIdx.x);

    curandState localState = state[i];

    if (i < N)
    {
        estimations[i] = (b - a) * f(a + (b - a) * curand_uniform(&localState));
        state[i] = localState;
    }
}

/**
 * Kernel para inicializar los estados de los números aleatorios a generar
 */
__global__ void Init(curandState *state)
{
    unsigned idy = threadIdx.x + (blockDim.x * blockIdx.x);
    curand_init(1234, idy, 0, &state[idy]);
}

int main(int argc, char **argv)
{
    char a_flag[] = "-a";
    char b_flag[] = "-b";
    char n_flag[] = "-n";

    double a, b;
    int n, BLOCK_SIZE = 1024;

    a = get_flag_value(argc, argv, a_flag);
    b = get_flag_value(argc, argv, b_flag);
    n = floor(get_flag_value(argc, argv, n_flag));

    double *cuda_random_estimations, *random_estimations;
    /**
     * Array con los estados para gen los números aleatorios
     */
    curandState *array_states;

    /**
     * inicializa los eventos
     */ 
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    /**
     * Reserva memoria para las estimaciones que se calculan en la GPU y que luego se copian en la CPU.
     */
    cudaMalloc((void **)&cuda_random_estimations, n * sizeof(double));
    cudaMalloc((void **)&array_states, n * sizeof(curandState));
    random_estimations = (double *)malloc(n * sizeof(double));

    printf("a = %f, b = %f, n = %d\n", a, b, n);

    double approx, real, error;

    int grid_size = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

    /**
     * Empieza a medir el tiempo que le toma al programa
     */
    cudaEventRecord(start);

    /**
     * Lanza el kernel para incializar cuRAND
     */
    Init<<<grid_size, BLOCK_SIZE>>>(array_states);

    /**
     * Lanza el kernel para calcular la integración Montecarlo
     */
    montecarlo_gpu<<<grid_size, BLOCK_SIZE>>>(cuda_random_estimations, array_states, a, b, n);

    /**
     * Copia los resultados calculados en la GPU (cuda_random_estimations) en la CPU (random_estimations)
     */
    cudaMemcpy(random_estimations, cuda_random_estimations, n * sizeof(double), cudaMemcpyDeviceToHost);
    /**
     * Sincroniza los kernels para que todos hayan acabado antes de continuar
     */
    cudaDeviceSynchronize();


    /**
     * Calcula el promedio de las estimaciones, el cual corresponde a la aproximación de la integral
     */
    for (int i = 0; i < n; i++)
    {
        approx += random_estimations[i];
    }

    approx /= n;

    /**
     * calcula la integral de la manera estándar
     */
    real = F(b) - F(a);
    /**
     * Obtiene el error de la aproximación
     */
    error = fabs((approx - real) / real);

    /**
     * Termina de medir el tiempo
     */
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    printf("Approximated value = %f (Tiempo empleado = %fs, %d iteraciones)\n", approx, milliseconds / 1000, n);
    printf("Real value = %f\n", real);
    printf("Error = %f\n", error);

    /**
     * Libera la memoria
     */
    cudaFree(array_states);
    cudaFree(cuda_random_estimations);
    free(random_estimations);

    return 0;
}