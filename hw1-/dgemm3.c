void dgemm3(double *C, double *A, double *B, int n) {
    // Use 16 registers for maximum register reuse
    double r[16];

    for (int i = 0; i < n; i += 4) {
        for (int j = 0; j < n; j += 4) {
            // Initialize the registers
            for (int k = 0; k < 16; k++) {
                r[k] = 0.0;
            }

            for (int k = 0; k < n; k++) {
                // Load data into registers
                register double a0 = A[i * n + k];
                register double a1 = A[(i + 1) * n + k];
                register double a2 = A[(i + 2) * n + k];
                register double a3 = A[(i + 3) * n + k];

                register double b0 = B[k * n + j];
                register double b1 = B[k * n + (j + 1)];
                register double b2 = B[k * n + (j + 2)];
                register double b3 = B[k * n + (j + 3)];

                // Perform the multiplications and accumulations
                r[0] += a0 * b0;
                r[1] += a0 * b1;
                r[2] += a0 * b2;
                r[3] += a0 * b3;

                r[4] += a1 * b0;
                r[5] += a1 * b1;
                r[6] += a1 * b2;
                r[7] += a1 * b3;

                r[8] += a2 * b0;
                r[9] += a2 * b1;
                r[10] += a2 * b2;
                r[11] += a2 * b3;

                r[12] += a3 * b0;
                r[13] += a3 * b1;
                r[14] += a3 * b2;
                r[15] += a3 * b3;
            }

            // Store the results back in C
            C[i * n + j] += r[0];
            C[i * n + j + 1] += r[1];
            C[i * n + j + 2] += r[2];
            C[i * n + j + 3] += r[3];

            C[(i + 1) * n + j] += r[4];
            C[(i + 1) * n + j + 1] += r[5];
            C[(i + 1) * n + j + 2] += r[6];
            C[(i + 1) * n + j + 3] += r[7];

            C[(i + 2) * n + j] += r[8];
            C[(i + 2) * n + j + 1] += r[9];
            C[(i + 2) * n + j + 2] += r[10];
            C[(i + 2) * n + j + 3] += r[11];

            C[(i + 3) * n + j] += r[12];
            C[(i + 3) * n + j + 1] += r[13];
            C[(i + 3) * n + j + 2] += r[14];
            C[(i + 3) * n + j + 3] += r[15];
        }
    }
}
