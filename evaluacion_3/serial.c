/**
 * @file serial.c
 * @author José Silva (mcvicksilvaa@gmail.com)
 * @brief Para ejecutar el archivo emplee el comando: gcc serial.c -o serial -lm && ./serial -a 0 -b 1 -n 100
 * @version 0.1
 * @date 2022-08-17
 *
 * @copyright Copyright (c) 2022
 *
 */

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
/**
 * Para medir el tiempo de ejecución
 */
#include <time.h>

/**
 * Consigue el valor de una flag proveniente de la linea de ocmandos
 *
 * @param argc cantidad de argumentos
 * @param argv array de argumentos
 * @param flag flag a parsear
 * @return int valor de la flag
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
double f(double x) {
    return pow(x, 2);
}

/**
 * Antiderivada de la función original
 */
double F(double x) {
    return pow(x, 3) / 3;
}

/**
 * Función para generar numeros aleatorios entre 0 y 1
 */
double rand_0_1() {
    return (double)rand() / (double)RAND_MAX;
}

/**
 * Función para calcular el método de integración Montecarlo
 */
double montecarlo_int(double a, double b, double n)
{
    double result = 0;
    for (int i = 0; i <= n; i++)
    {
        result += (b - a) * f(a + (b - a) * rand_0_1());
    }
    result /= n;
    return result;
}

int main(int argc, char **argv)
{
    /**
     * Garantiza una semilla distinta para ejecución del programa
     */
    srand(time(NULL));
    char a_flag[] = "-a";
    char b_flag[] = "-b";
    char n_flag[] = "-n";

    double a, b, n;
    double X;

    a = get_flag_value(argc, argv, a_flag);
    b = get_flag_value(argc, argv, b_flag);
    n = get_flag_value(argc, argv, n_flag);

    printf("a = %f, b = %f, n = %f\n", a, b, n);

    double approx, real, error;

    /**
     * Mide el tiempo de ejecución de la aproximación
     */
    clock_t begin = clock();
    /**
     * Calcula la integral con Montecarlo
     */
    approx = montecarlo_int(a, b, n);
    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

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