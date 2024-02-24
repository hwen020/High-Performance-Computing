#ifndef __MY_C__
#define __MY_C__

#include "include.h"

int mydgetrf(double *A, int *ipiv, int n) {
    // TODO
    // The return value (an integer) can be 0 or 1
    // If 0, the matrix is irreducible, and the result will be ignored
    // If 1, the result is valid

    /* add your code here */
    int i, j, k;
    int max_Aii;
    double max;
    double* temp_Aii = (double*)malloc(n * sizeof(double)); // a temp row for swap
    
    for (i = 0; i < n; i++) {
        // find the biggest |A[i,i]|
        max = fabs(A[i*n + i]);
        max_Aii = i; 
        for (j = i+1; j < n; j++) {
            if ( fabs(A[j*n + i]) > max ) { // update the biggest element in i-th col
                max = fabs(A[j*n + i]);
                max_Aii = j; // |A[j,i]| is the biggest in rest of i-th col
            }
        }
        // swap row i and max_Aii and the corrosponding vector
        if (max == 0) {
            // printf("WARNING! A is singular, or near so... \n");
            return -1;
        } else {
            if (max_Aii != i) { 
                int temp = ipiv[i];
                ipiv[i] = ipiv[max_Aii];
                ipiv[max_Aii] = temp;
                memcpy(temp_Aii, A + i*n, n*sizeof(double));
                memcpy(A + i*n, A + max_Aii*n, n*sizeof(double));
                memcpy(A + max_Aii*n, temp_Aii, n*sizeof(double));
            }
        }

        for (j = i+1; j < n; j++) {
            A[j*n + i] = A[j*n+i] / A[i*n + i];
            for (k = i+1; k < n; k++) {
                A[j*n + k] -= A[j*n + i] * A[i*n + k];
            }
        }
    }
    free(temp_Aii);

    return 1; // Successful LU factorization.
}

void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv) {
    // // TODO
    int i, j;
    double* y = (double*)malloc(n * sizeof(double)); // vector y
    
    if (UPLO == 'L') { // lower triangle, calculating vector y by L*y=b
        y[0] = B[ipiv[0]]; // y0*1 + y1*0 + y2*0+...+y(n-1)*0=b0
        for (i = 1;i < n;i++) {
            double sum = 0.0;
            for (j = 0; j < i; j++) {  // only use the lower triangle of A
                sum += A[i*n + j] * y[j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    } 
    if (UPLO == 'U') { 
        y[n-1] = B[n-1]/A[n*n-1];// y[n-1] is correct!!!
        for (i = n-2;i >= 0;i--) { // the other y[i] is incorrect!!!
            double sum = 0.0;
            for (j = i+1;j < n;j++) {
                sum += A[i*n + j] * y[j];
            }
            y[i] = (B[i] - sum)/A[i*n + i];
        }
    }
    memcpy(B, y, n*sizeof(double)); // copy back to vector B
    free(y);
    return;
}

void my_f(double *A, double *B, int n) {
    int *ipiv = (int *)malloc(n * sizeof(int));
    if (ipiv == NULL) {
        printf("Memory allocation failed for ipiv.\n");
        return;
    }

    if (mydgetrf(A, ipiv, n) == 0) {
        printf("LU factorization failed: coefficient matrix is singular.\n");
        free(ipiv);
        return;
    }

    mydtrsv('L', A, B, n, ipiv);
    mydtrsv('U', A, B, n, ipiv);

    free(ipiv);
}

#endif
