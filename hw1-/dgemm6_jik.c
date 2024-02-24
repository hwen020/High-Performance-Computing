void dgemm6_jik(double *C,double *A,double *B,int n)
{
    int j, i, k;
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
        {
            register double sum = 0.0;
            for (k = 0; k < n; k++)
                sum += A[i * n + k] * B[k * n + j];
            C[i * n + j] = sum;
        }
}

void dgemm6_jik2(double *C,double *A,double *B,int n)
{
    int j, jj, i, ii, k, kk;
    int b = 1; //change b to the number you want
    for (j = 0; j < n; j += b)
        for (i = 0; i < n; i += b)
            for (k = 0; k < n; k += b)
                for (jj = j; jj < j + b; jj++)
                    for (ii = i; ii < i + b; ii++)
                    {
                        register double sum = 0.0;
                        for (kk = k; kk < k + b; kk++)
                            sum += A[ii * n + kk] * B[kk * n + jj];
                        C[ii * n + jj] = sum;
                    }
}