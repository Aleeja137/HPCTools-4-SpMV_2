#include <stdio.h>
#include <stdlib.h>
#include <mkl_spblas.h>

int main() {
    // Matrix and vector size (assuming square matrix)
    int size = 4;  // Replace with your matrix size

    // Allocate and initialize dense matrix and vector
    double *mat = (double *)malloc(size * size * sizeof(double));
    double *vec = (double *)malloc(size * sizeof(double));
    double *result = (double *)malloc(size * sizeof(double));
    
    // Fill in the matrix and vector with sample values
    // Replace this with actual data
    for (int i = 0; i < size * size; i++) mat[i] = (i % 3 == 0) ? i + 1.0 : 0.0; // sparse pattern
    for (int i = 0; i < size; i++) vec[i] = 1.0;

    // Create MKL sparse matrix in COO format
    sparse_matrix_t A;
    int nnz = 0;  // Count non-zeros
    for (int i = 0; i < size * size; i++) if (mat[i] != 0.0) nnz++;

    // Allocate COO format arrays
    MKL_INT *rows = (MKL_INT *)malloc(nnz * sizeof(MKL_INT));
    MKL_INT *cols = (MKL_INT *)malloc(nnz * sizeof(MKL_INT));
    double *values = (double *)malloc(nnz * sizeof(double));

    // Fill COO format arrays
    int idx = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (mat[i * size + j] != 0.0) {
                rows[idx] = i;
                cols[idx] = j;
                values[idx] = mat[i * size + j];
                idx++;
            }
        }
    }

    // Create sparse matrix handle in MKL COO format
    sparse_status_t status = mkl_sparse_d_create_coo(&A, SPARSE_INDEX_BASE_ZERO, size, size, nnz, rows, cols, values);
    if (status != SPARSE_STATUS_SUCCESS) {
        printf("Error creating sparse matrix\n");
        return -1;
    }

    // Define matrix descriptor
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Perform sparse matrix-vector multiplication
    double alpha = 1.0;
    double beta = 0.0;
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, vec, beta, result);
    if (status != SPARSE_STATUS_SUCCESS) {
        printf("Error performing sparse matrix-vector multiplication\n");
        return -1;
    }

    // Output result
    printf("Result vector:\n");
    for (int i = 0; i < size; i++) {
        printf("%f\n", result[i]);
    }

    // Cleanup
    mkl_sparse_destroy(A);
    free(mat);
    free(vec);
    free(result);
    free(rows);
    free(cols);
    free(values);

    return 0;
}
