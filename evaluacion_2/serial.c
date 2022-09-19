/**
 * @file serial.c
 * @author José Silva (mcvicksilvaa@gmail.com)
 * @brief Para ejecutar el archivo emplee el comando: gcc serial.c -o serial -lm && ./serial -nEXP 5000 -nPI 5000 -x 1
 * @version 0.1
 * @date 2022-08-17
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <stdio.h>
#include <stdlib.h>
/**
 * Para utilizar las funciones exp y pow
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
 * Valor de π
 */
double PI = 3.14159265358979323846;

/**
 * Consigue el valor de una flag proveniente de la linea de ocmandos
 * 
 * @param argc cantidad de argumentos
 * @param argv array de argumentos
 * @param flag flag a parsear
 * @return int valor de la flag
 */
int get_flag_value(int argc, char **argv, char flag[]) {
    int n = -1, arg = 0;
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
 * @brief Aproxima el valor de π utilizando la formula de Leibniz
 *
 * @param n cantidad de términos
 * @return double
 */
double compute_pi(int n) {
    double pi = 0;
    /**
     * valor alternante
     */
    double alternate = 1;
    for (int i = 0; i < n; i++)
    {
        /**
         * alterna entre los inversos de los números impares
         */
        pi += alternate / (2 * i + 1);
        alternate *= -1;
    }
    return 4 * pi;
}

/**
 * @brief Calcula el factorial de un número recursivamente
 * 
 * @param n número al cual se le va a calcular el factorial
 * @return double 
 */
double fact(double n) {
    if (n >= 1) {
        return n * fact(n - 1);
    } else {
        return 1;
    }
}

/**
 * @brief Esencialmente calcula lo mismo que exp(x), pero mediante serie de potencias
 *
 * @param x valor de x que se va a evaluar
 * @param n_exp cantidad de términos para la seri de potencias
 * @param n_pi cantidad de términos para aproximar π
 * @return double
 */
double exponential(double x, int n_exp, int n_pi) {
    double expr = 0;
    double pi;
    pi = compute_pi(n_pi);

    for (int k = 0; k < n_exp; k++)
    {
        expr += pow(x, k) / (fact(k));
    }

    return pow(expr, pi);
}

int main(int argc, char **argv) {
    char exp_n_terms_flag[] = "-nEXP";
    char pi_n_terms_flag[] = "-nPI";
    char x_flag[] = "-x";

    int n_exp, n_pi;
    double X;

    n_exp = get_flag_value(argc, argv, exp_n_terms_flag);
    n_pi = get_flag_value(argc, argv, pi_n_terms_flag);
    if (strcmp(x_flag, argv[argc - 2]) != 0) {
        printf("Missing '-x' flag.\n");
        return 0;
    }

    X = atof(argv[argc - 1]);
 
    double approx, real, error;

    /**
     * Mide el tiempo de ejecución de la aproximación
     */
    clock_t begin = clock();
    approx = exponential(X, n_exp, n_pi);
    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    /**
     * calcula e^(πx) con funciones estándares
     */
    real = pow(exp(X), PI);
    /**
     * Obtiene el error de la aproximación
     */
    error = fabs((approx - real) / real);

    printf("Approximated value = %f (Tiempo empleado = %fs)\n", approx, time_spent);
    printf("Real value = %f\n", real);
    printf("Error = %f\n", error);

    return 0;
}