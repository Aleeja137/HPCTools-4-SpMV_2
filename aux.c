#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spblas.h>
#include <mkl.h>                // Intel MKL for CBLAS and Sparse BLAS
#include <mkl_spblas.h>
#include "timer.h"              // Assumed user-defined timing library
#include "spmv.h"               // Assumed user-defined header for custom matrix functions

#define DEFAULT_SIZE 4
#define DEFAULT_DENSITY 1

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed) {
    unsigned int nnz = 0;
    srand(seed);
    for (unsigned int i = 0; i < n * n; i++) {
        if ((rand() % 100) / 100.0 < density) {
            mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
            nnz++;
        } else {
            mat[i] = 0;
        }
    }
    return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed) {
    srand(seed);
    for (unsigned int i = 0; i < size; i++) {
        vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }
    return 0;
}

// Sparse matrix-vector multiplication using MKL (CSR format)
void sparse_mkl_matvec(double *values, int *columns, int *rowIndex, double *x, double *y, int n) {
    sparse_matrix_t mkl_A;
    struct matrix_descr descr;

    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_d_create_csr(&mkl_A, SPARSE_INDEX_BASE_ZERO, n, n, rowIndex, rowIndex + 1, columns, values);
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_A, descr, x, 0.0, y);
    mkl_sparse_destroy(mkl_A);
}

// Main function
int main(int argc, char *argv[]) {
    int n = DEFAULT_SIZE;
    double density = DEFAULT_DENSITY;
    unsigned int seed = 0;

    double *dense_matrix = (double *)malloc(n * n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    double *y_gnu = (double *)malloc(n * sizeof(double));
    double *y_mkl = (double *)malloc(n * sizeof(double));

    // Populate matrices and vectors
    populate_sparse_matrix(dense_matrix, n, density, seed);
    populate_vector(x, n, seed);

    // Timing GNU CBLAS dense matrix-vector multiplication
    timer_start();
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, dense_matrix, n, x, 1, 0.0, y_gnu, 1);
    timer_stop();
    printf("GNU CBLAS dense matrix-vector multiplication time: %f ms\n", timer_elapsed());

    // Timing MKL dense matrix-vector multiplication
    timer_start();
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, y, 1);
    timer_stop();
    printf("MKL dense matrix-vector multiplication time: %f ms\n", timer_elapsed());

    // Convert dense matrix to CSR for sparse matrix operations
    int nnz = 0;    // Count of non-zero elements
    int *rowIndex = (int *)malloc((n + 1) * sizeof(int));
    int *columns = (int *)malloc(n * n * sizeof(int));
    double *values = (double *)malloc(n * n * sizeof(double));

    rowIndex[0] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dense_matrix[i * n + j] != 0) {
                values[nnz] = dense_matrix[i * n + j];
                columns[nnz] = j;
                nnz++;
            }
        }
        rowIndex[i + 1] = nnz;
    }

    // Timing MKL sparse matrix-vector multiplication
    timer_start();
    sparse_mkl_matvec(values, columns, rowIndex, x, y_mkl, n);
    timer_stop();
    printf("MKL sparse matrix-vector multiplication time: %f ms\n", timer_elapsed());

    // Clean up memory
    free(dense_matrix);
    free(x);
    free(y_gnu);
    free(y_mkl);
    free(rowIndex);
    free(columns);
    free(values);

    return 0;
}
