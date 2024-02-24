void dgemm7(double *C, double *A, double *B, int n) {
    // Experiment with different block sizes for cache and register blocking
    int cacheBlockSize = 64; // Try 64, 128 first. Range of 16 to 256
    int registerBlockSize = 4; // Try 4, 8, 16

    for (int i0 = 0; i0 < n; i0 += cacheBlockSize) {
        for (int j0 = 0; j0 < n; j0 += cacheBlockSize) {
            for (int k0 = 0; k0 < n; k0 += cacheBlockSize) {
                for (int i = i0; i < i0 + cacheBlockSize; i += registerBlockSize) {
                    for (int j = j0; j < j0 + cacheBlockSize; j += registerBlockSize) {
                        for (int k = k0; k < k0 + cacheBlockSize; k += registerBlockSize) {
                            // Initialize registers for the result block
                            double resultRegisters[registerBlockSize][registerBlockSize];
                            for (int ii = 0; ii < registerBlockSize; ii++) {
                                for (int jj = 0; jj < registerBlockSize; jj++) {
                                    resultRegisters[ii][jj] = 0.0;
                                }
                            }

                            // Perform matrix multiplication within registers
                            for (int kk = 0; kk < registerBlockSize; kk++) {
                                for (int i_inner = 0; i_inner < registerBlockSize; i_inner++) {
                                    for (int j_inner = 0; j_inner < registerBlockSize; j_inner++) {
                                        resultRegisters[i_inner][j_inner] += A[(i + i_inner) * n + (k + kk)] *
                                            B[(k + kk) * n + (j + j_inner)];
                                    }
                                }
                            }

                            // Accumulate the result in the main C matrix
                            for (int i_inner = 0; i_inner < registerBlockSize; i_inner++) {
                                for (int j_inner = 0; j_inner < registerBlockSize; j_inner++) {
                                    C[(i + i_inner) * n + (j + j_inner)] += resultRegisters[i_inner][j_inner];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
