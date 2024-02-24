#ifndef __SIEVE2_C__
#define __SIEVE2_C__

#include "include.h"

void sieve2(unsigned long long *global_count,unsigned long long n,int pnum,int pid)
{
    unsigned long long low_value = 3 + pid * ((n - 1) / pnum); // Start from 3 to skip even numbers
    unsigned long long high_value = 1 + (pid + 1) * ((n - 1) / pnum);
    unsigned long long size = (high_value - low_value + 1) / 2; // Divide size by 2 for odd numbers only

    if (1 + (n - 1) / pnum < (int)sqrt((double)n))
    {
        if (pid == 0)
            printf("Error: Too many processes.\n");
        MPI_Finalize();
        exit(0);
    }

    char *marked = (char *)malloc(size); // Array for marking multiples. 1 means multiple and 0 means prime
    if (marked == NULL)
    {
        printf("Error: Cannot allocate enough memory.\n");
        MPI_Finalize();
        exit(0);
    }
    memset(marked, 0, size);

    unsigned long long index = 0; // Index of the current prime among all primes (only works for process 0)
    unsigned long long prime = 3; // Start with the first odd prime
    unsigned long long first; // Index of the first multiple among values handled by this process

    do
    {
        if (prime * prime > low_value)
            first = (prime * prime - low_value) / 2; // Divide by 2 for odd numbers only
        else
            if (low_value % prime == 0)
                first = 0;
            else
                first = (prime - low_value % prime) % 2 == 0 ? (prime - low_value % prime) / 2 : (prime - low_value % prime + prime) / 2;

        for (unsigned long long i = first; i < size; i += prime)
            marked[i] = 1;
        index++;

        while (marked[index] == 1)
            index++;
        prime = 2 * index + 3; // Adjust for odd primes only
    } while (prime * prime <= n);

    unsigned long long count = 0; // Local count of primes
    if (pid == 0)
        count++; // 2 is prime
    for (unsigned long long i = 0; i < size; i++)
        if (marked[i] == 0)
            count++;

    MPI_Reduce(&count, global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    free(marked);
}

#endif