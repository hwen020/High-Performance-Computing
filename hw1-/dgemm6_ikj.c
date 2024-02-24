void dgemm6_ikj(double *C,double *A,double *B,int n)
{
    int i, k, j;
    for (i = 0; i < n; i++)
        for (k = 0; k < n; k++)
        {
            register double r = A[i * n + k];
            for (j = 0; j < n; j++)
                C[i * n + j] += r * B[k * n + j];
        }
}

void dgemm6_ikj2(double *C,double *A,double *B,int n)
{
    int i, ii, k, kk, j, jj;
    int b = 1; //change b to the number you want
    for (i = 0; i < n; i += b)
        for (k = 0; k < n; k += b)
            for (j = 0; j < n; j += b)
                for (ii = i; ii < i + b; ii++)
                    for (kk = k; kk < k + b; kk++)
                    {
                        register double r = A[ii * n + kk];
                        for (jj = j; jj < j + b; jj++)
                            C[ii * n + jj] += r * B[kk * n + jj];
                    }
}