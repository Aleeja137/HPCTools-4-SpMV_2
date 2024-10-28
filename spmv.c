#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spblas.h>
#include "timer.h"
#include "spmv.h"
#include "mkl_spblas.h"

// #define DEFAULT_SIZE 16384
// #define DEFAULT_DENSITY 0.1

#define DEFAULT_SIZE 4
#define DEFAULT_DENSITY 1

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
  unsigned int nnz = 0;

  srand(seed);

  for (unsigned int i = 0; i < n * n; i++) {
    if ((rand() % 100) / 100.0 < density) {
      // Get a pseudorandom value between -9.99 e 9.99
      mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
      nnz++;
    } else {
      mat[i] = 0;
    }
  }

  return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
  srand(seed);

  for (unsigned int i = 0; i < size; i++) {
    vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
  }

  return size;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return fabs(x - y) <= epsilon * fabs(x);
  // see Knuth section 4.2.2 pages 217-218
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    if (!is_nearly_equal(ref[i], result[i]))
      return 0;
  }

  return 1;
}

unsigned int check_result_gsl(gsl_vector* ref, gsl_vector* result, unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    double ref_elem = gsl_vector_get(ref, i);
    double result_elem = gsl_vector_get(result, i);
    if (!is_nearly_equal(ref_elem, result_elem))
      return 0;
  }

  return 1;
}


int main(int argc, char *argv[])
{
  int size, i, j;  // number of rows and cols (size x size matrix)
  double density;  // aprox. ratio of non-zero values

  if (argc < 2) {
    size = DEFAULT_SIZE;
    density = DEFAULT_DENSITY;
  } else if (argc < 3) {
    size = atoi(argv[1]);
    density = DEFAULT_DENSITY;
  } else {
    size = atoi(argv[1]);
    density = atoi(argv[2]);
  }

  double *mat, *vec, *refsol, *mysol;

  // Reserve memory
  mat = (double *) malloc(size * size * sizeof(double));
  vec = (double *) malloc(size * sizeof(double));
  refsol = (double *) malloc(size * sizeof(double));
  mysol = (double *) malloc(size * sizeof(double));

  // Initialize matrix and vector
  unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
  populate_vector(vec, size, 2);

  // Print initial control information
  printf("Matriz size: %d x %d (%d elements)\n", size, size, size*size);
  printf("%d non-zero elements (%.2lf%%)\n", nnz, (double) nnz / (size*size) * 100.0);

  #ifdef GCC
    printf("Compiled with gcc compiler\n\n");
  #elif ICC
      printf("Compiled with icc compiler\n\n");
  #endif


  /* --- Dense computation using CBLAS (eg. GSL's CBLAS implementation) --- */
  #ifdef GCC
    printf("Dense computation\n----------------\n");
    timeinfo start, now;
    timestamp(&start);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

    timestamp(&now);
    printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));
  #elif ICC
    printf("Dense computation\n----------------\n");
    timeinfo start, now;
    timestamp(&start);

    // TODO - change for MKL's blas multiplication
    cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

    timestamp(&now);
    printf("Time taken by MLK dense computation: %ld ms\n", diff_milli(&start, &now));
  #endif


  /* --- Dense computation using your own dense implementation --- */
  timestamp(&start);

  my_dense(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  // Check results
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

 
  /* --- Sparse Matrix - Dense Vector computation --- */

  // First initialize like a 'dense' matrix (nzmax = size*size) with COO format
  gsl_vector *gsl_vec = gsl_vector_alloc(size);
  gsl_vector *gsl_vec_result_coo = gsl_vector_alloc(size);
  gsl_vector *gsl_vec_result_csr = gsl_vector_alloc(size);
  gsl_vector *gsl_vec_result_csc = gsl_vector_alloc(size);
  gsl_vector *gsl_vec_result_mine = gsl_vector_calloc(size);
  gsl_spmatrix *smat_coo = gsl_spmatrix_alloc_nzmax(size, size, size*size, GSL_SPMATRIX_COO);

  // Set sparse matrix values
  for (i=0; i<size; i++)
  {
    gsl_vector_set(gsl_vec, i, vec[i]);
    for(j=0; j<size; j++)
      if (mat[i*size+j]!=0)
        gsl_spmatrix_set(smat_coo, i, j, mat[i*size+j]);
  }  

  // Compress into different formats
  smat_coo = gsl_spmatrix_compress(smat_coo, GSL_SPMATRIX_COO);
  gsl_spmatrix *smat_csr = gsl_spmatrix_compress(smat_coo, GSL_SPMATRIX_CSR);
  gsl_spmatrix *smat_csc = gsl_spmatrix_compress(smat_coo, GSL_SPMATRIX_CSC);

  /* --- Sparse computation using GSL's sparse algebra functions --- */
  #ifdef GCC
    printf("\nSparse computation\n----------------\n");
    timestamp(&start);
    gsl_spblas_dgemv(CblasNoTrans, 1.0, smat_coo, gsl_vec, 0.0, gsl_vec_result_coo);
    timestamp(&now);
    printf("Time taken by GSL COO sparse computation: %ld ms\n", diff_milli(&start, &now));

    timestamp(&start);
    gsl_spblas_dgemv(CblasNoTrans, 1.0, smat_csr, gsl_vec, 0.0, gsl_vec_result_csr);
    timestamp(&now);
    printf("Time taken by GSL CSR sparse computation: %ld ms\n", diff_milli(&start, &now));

    timestamp(&start);
    gsl_spblas_dgemv(CblasNoTrans, 1.0, smat_csc, gsl_vec, 0.0, gsl_vec_result_csc);
    timestamp(&now);
    printf("Time taken by GSL CSC sparse computation: %ld ms\n", diff_milli(&start, &now));

  /* --- Sparse computation using MKL's sparse algebra functions --- */
  #elif ICC
    // TODO - implement the same using mkl_sparse_d_mv function. 
  #endif

  /* --- Your own sparse implementation --- */

  // COO format computation
  timestamp(&start);
  my_sparse_coo(smat_coo, gsl_vec, gsl_vec_result_mine);
  timestamp(&now);
  if (check_result_gsl(gsl_vec_result_coo, gsl_vec_result_mine, size) == 1)
    printf("COO result is ok! - ");
  else
    printf("COO result is wrong! - ");
  printf("Time taken by my COO sparse matrix-vector product: %ld ms\n", diff_milli(&start, &now));
  gsl_vec_result_mine = gsl_vector_calloc(size);

  // CSR format computation
  timestamp(&start);
  my_sparse_csr(size, smat_csr, gsl_vec, gsl_vec_result_mine);
  timestamp(&now);
  if (check_result_gsl(gsl_vec_result_csr, gsl_vec_result_mine, size) == 1)
    printf("CSR result is ok! - ");
  else
    printf("CSR result is wrong! - ");
  printf("Time taken by my CSR sparse matrix-vector product: %ld ms\n", diff_milli(&start, &now));
  gsl_vec_result_mine = gsl_vector_calloc(size);

  // CSC format computation
  timestamp(&start);
  my_sparse_csc(size, smat_csc, gsl_vec, gsl_vec_result_mine);
  timestamp(&now);
  if (check_result_gsl(gsl_vec_result_csc, gsl_vec_result_mine, size) == 1)
    printf("CSC result is ok! - ");
  else
    printf("CSC result is wrong! - ");
  printf("Time taken by my CSC sparse matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  
  // Free resources
  free(mat);
  free(vec);
  free(refsol);
  free(mysol);
  gsl_spmatrix_free(smat_csr);
  gsl_spmatrix_free(smat_coo);
  gsl_vector_free(gsl_vec);
  gsl_vector_free(gsl_vec_result_coo);
  gsl_vector_free(gsl_vec_result_csr);
  gsl_vector_free(gsl_vec_result_csc);
  gsl_vector_free(gsl_vec_result_mine);

  return 0;
}
