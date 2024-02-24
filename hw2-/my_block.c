#ifndef __MY_BLOCK_C__
#define __MY_BLOCK_C__

#include "include.h"
#include <string.h>

void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b) {
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    int i1, j1, k1;

    for (i1 = i; i1 < i+b; i1 += 2) {
        for (k1 = k; k1 < k+b; k1 += 2) {
            // reg for A
            register int ra1 = i1*n + k1;// A[i1, k1]
            register int ra2 = ra1 + n;// A[(i1+1), k1]
            register double a00 = A[ra1];// A[i1, k1]
            register double a01 = A[ra1+1];// A[i1, k1+1]
            register double a10 = A[ra2];// A[(i1+1), k1]
            register double a11 = A[ra2+1];// A[(i1+1), (k1+1)]
            for (j1 = j; j1 < j+b; j1 += 2) {
                // reg for C
                register int rc1 = i1*n + j1;// C[i1, j1]
                register int rc2 = rc1 + n;// C[i1+1, j1]
                register double c00 = C[rc1];// C[i1, j1]
                register double c01 = C[rc1+1];// C[i1, j1+1]
                register double c10 = C[rc2];// C[i1+1, j1]
                register double c11 = C[rc2+1];// C[i1+1, j1+1]
                // reg for B
                register int rb1 = k1*n + j1;// B[k1, j1]
                register int rb2 = rb1 + n;// B[(k1+1), j1]
                register double b00 = B[rb1];// B[k1, j1]
                register double b01 = B[rb1+1];// B[k1, j1+1]
                register double b10 = B[rb2];// B[(k1+1), j1]
                register double b11 = B[rb2+1];// B[k1+1, j1+1]
                // calculate C
                c00 -= a00 * b00 + a01 * b10;
                c01 -= a00 * b01 + a01 * b11;
                c10 -= a10 * b00 + a11 * b10;
                c11 -= a10 * b01 + a11 * b11;
                C[rc1] = c00;
                C[rc1+1] = c01;
                C[rc2] = c10;
                C[rc2+1] = c11;
            }
        }
    }
}


int mydgetrf_block(double *A, int *ipiv, int n) {
    int b=4; //block size
    int i, j, k;
    int max_Aii;
    int ib;
    double max;
    double* temp_Aii = (double*)malloc(n * sizeof(double)); // a temp row for swap
    
    for (ib = 0;ib < n;ib += b) {
        for (i = ib; i < ib+b; i++) {
            max = fabs(A[i*n + i]);
            max_Aii = i;
            for (j = i+1;j < n;j++) {
                if ( fabs(A[j*n + i]) > max ) {
                    max = fabs(A[j*n + i]);
                    max_Aii = j; // |A[j,i]| is the biggest in rest of i-th col
                }
            }
            // swap row i and max_Aii and the corrosponding vector
            if (max == 0) {
                return -1;
            } else {
                if (max_Aii != i) { // bug here!!!!
                    // swap i-th and max_Aii-th element of vector pivot
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[max_Aii];
                    ipiv[max_Aii] = temp;
                    // swap i-th row and max_Aii-th row
                    memcpy(temp_Aii, A + i*n, n*sizeof(double));
                    memcpy(A + i*n, A + max_Aii*n, n*sizeof(double));
                    memcpy(A + max_Aii*n, temp_Aii, n*sizeof(double));
                }
            }
            for (j = i+1; j < n; j++) {
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                for (k = i+1; k < ib+b; k++) {
                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }
            }
        }
        // A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n), update next b rows of U
        for (i = ib; i < ib+b; i++) {
            for (j = ib+b;j < n;j++) { // j:[end+1=ib+b-1+1, n]
                double sum = 0.0;
                for (k = ib;k < i;k++) {
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }

        // update A(end+1=ib+b-1+1 : n, end+1:n) block by block
        for (i = ib+b;i < n;i+=b) {
            for (j = ib+b;j < n;j+=b) {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }

    return 1;
}


void my_block_f(double *A, double *B, int n) {
    int *ipiv = (int *)malloc(n * sizeof(int));
    if (ipiv == NULL) {
        printf("Memory allocation failed for ipiv.\n");
        return;
    }

    if (mydgetrf_block(A, ipiv, n) == 0) {
        printf("LU factorization failed: coefficient matrix is singular.\n");
        free(ipiv);
        return;
    }

    mydtrsv('L', A, B, n, ipiv);
    mydtrsv('U', A, B, n, ipiv);

    free(ipiv);
}

#endif
