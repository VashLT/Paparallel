/**
 * @file prime_multi.c
 * @author José Silva (mcvicksilvaa@gmail.com)
 * @brief Para ejecutar el archivo emplee el comando: mpicc -lm ./prime_multi.c; mpirun -np 2 a.out 100
 * @version 0.1
 * @date 2022-07-10
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/**
 * @brief Algoritmo para determinar si un número dado es primo o no
 * 
 * @param number 
 * @return int 
 */
int is_prime(int number)
{
    if (number == 2 || number == 3)
    {
        return 1;
    }

    if (number <= 1 || number % 2 == 0 || number % 3 == 0)
    {
        return 0;
   }

    for (int i = 5; i * i <= number; i += 6)
    {
        if (number % i == 0 || number % (i + 2) == 0)
            return 0;
    }

    return 1;
}

int main(int argc, char **argv)
{
    int number;

    number = atoi(argv[1]);

    if (number < 2)
    {
        printf("The input number must be greater or equal than 2.");
        return 0;
    }

    /**
     * curr_rank: rank actual, se le envian procesos de computo
     * prime: número a determinar si es primo o no
     * index: index para almacenar el array de primos
     */
    int rank, size, prime, num_is_prime, curr_rank = 0, index = 0;
    int primes[number];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /**
     * se asume que se necesitan minimo 2 ranks para ejecutar el algoritmo
     */
    if (size < 2)
    {
        fprintf(stderr, "World size must be greater than 1 for %s\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        for (int i = 2; i < number; i += 1)
        {
            curr_rank++;
            if (curr_rank >= size) {
                curr_rank = 1;
            }
            MPI_Send(
                &i,
                1,
                MPI_INT,
                curr_rank,
                0,
                MPI_COMM_WORLD);

            MPI_Recv(
                &num_is_prime,
                1,
                MPI_INT,
                curr_rank,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
                
            if (num_is_prime) {
                primes[index] = i;
                index++;
            }
        }

        for (int j = 0; j < index; j += 1) {
            printf("%d is Prime\n", primes[j]);
        }
    } else {
        int candidate;
        for (int i = 2 + (rank - 1); i < number; i += (size - 1))
        {
            MPI_Recv(
                &candidate,
                1,
                MPI_INT,
                0,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

            num_is_prime = is_prime(candidate);

            MPI_Send(
                &num_is_prime,
                1,
                MPI_INT,
                0,
                0,
                MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();

    return 0;
}