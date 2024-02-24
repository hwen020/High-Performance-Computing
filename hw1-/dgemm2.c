void dgemm2(double *C, double *A, double *B, int n) {
    // Allocate 12 registers for a, b, and c
    register double c00, c01, c10, c11, c20, c21, c30, c31, a00, a10, a01, a11, b00, b01, b10, b11;

    for (int i = 0; i < n; i += 2) {
        for (int j = 0; j < n; j += 2) {
            int t = i * n + j;
            int tt = t + n;
            c00 = C[t];
            c01 = C[t + 1];
            c10 = C[tt];
            c11 = C[tt + 1];

            int t2 = (i + 1) * n + j;
            int t2t = t2 + n;
            c20 = C[t2];
            c21 = C[t2 + 1];
            c30 = C[t2t];
            c31 = C[t2t + 1];

            for (int k = 0; k < n; k += 2) {
                int ta = i * n + k;
                int tta = ta + n;
                int ta2 = (i + 1) * n + k;
                int tta2 = ta2 + n;  
                int tb = k * n + j;
                int ttb = tb + n;
                int tb2 = k * n + j;
                int tb2t = tb2 + n;

                a00 = A[ta];
                a10 = A[tta];
                a01 = A[ta2];
                a11 = A[tta2];  
                b00 = B[tb];
                b01 = B[tb + 1];

                // Perform the 2x2 matrix multiplication
                c00 += a00 * b00;
                c01 += a01 * b00;
                c10 += a00 * b01;
                c11 += a01 * b01;

                c20 += a10 * b00;
                c21 += a11 * b00;
                c30 += a10 * b01;
                c31 += a11 * b01;

                a00 = A[tta + 1];
                a10 = A[tta2 + 1];  
                a01 = A[tta + 1];
                a11 = A[tta2 + 1];  
                b00 = B[tb2t];
                b01 = B[tb2t + 1];

                // Continue with the multiplication
                c00 += a00 * b00;
                c01 += a01 * b00;
                c10 += a00 * b01;
                c11 += a01 * b01;

                c20 += a10 * b00;
                c21 += a11 * b00;
                c30 += a10 * b01;
                c31 += a11 * b01;
            }

            // Store the results back in C
            C[t] = c00;
            C[t + 1] = c01;
            C[tt] = c10;
            C[tt + 1] = c11;
            C[t2] = c20;
            C[t2 + 1] = c21;
            C[t2t] = c30;
            C[t2t + 1] = c31;
        }
    }
}
