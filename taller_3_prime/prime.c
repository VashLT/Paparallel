#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

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

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // printf("Rank %d:\n", rank + 1);
    // 2 is the first prime number
    for (int i = 2 + rank; i < number; i+=size) {
        // printf("    test: %d\n", i);
        if (is_prime(i))
        {
            printf("    %d\n", i);
        }
    }
    MPI_Finalize();

    return 0;
}