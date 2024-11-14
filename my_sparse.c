#include "spmv.h"
#include <stddef.h>
#ifdef GCC
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>


void my_sparse_coo(gsl_spmatrix *mat, gsl_vector *vec, gsl_vector *result)
{
  size_t i, tot_elem = mat->nz;
  int row, col;

  for(i=0; i<tot_elem; i++)
  {
    row = mat->i[i];
    col = mat->p[i];
    // Do the multiplication and cumulative sum of result on target position
    double vec_element = gsl_vector_get(vec, col);
    double mat_element = gsl_spmatrix_get(mat, row, col);
    double result_element = gsl_vector_get(result, row);
    
    // This works since gsl_vector_calloc initializes vector to zero, otherwise an initial set whould be done
    result_element += mat_element * vec_element;
    gsl_vector_set(result, row, result_element);
  }
}

void my_sparse_csr(size_t n, gsl_spmatrix *mat, gsl_vector *vec, gsl_vector *result)
{
  size_t i;
  int j, row, col;
  int n_elem, n_elem_total = 0;

  for(i=0; i<n; i++)
  {
    // Get how many elements are in the row i
    n_elem = mat->p[i+1] - mat->p[i];

    // For each element, get their column index
    for (j=0; j<n_elem; j++)
    {
      // The only difference with CSC lies here
      row = i;
      col = mat->i[n_elem_total+j];

      // Do the multiplication and cumulative sum of result on target position
      double vec_element = gsl_vector_get(vec, col);
      double mat_element = gsl_spmatrix_get(mat, row, col);
      double result_element = gsl_vector_get(result, row);

      // This works since gsl_vector_calloc initializes vector to zero, otherwise an initial set whould be done
      result_element += mat_element * vec_element;
      gsl_vector_set(result, row, result_element);
    }
    n_elem_total += n_elem;
  }
}

void my_sparse_csc(size_t n, gsl_spmatrix *mat, gsl_vector *vec, gsl_vector *result)
{
  size_t i;
  int j, row, col;
  int n_elem, n_elem_total = 0;

  for(i=0; i<n; i++)
  {
    // Get how many elements are in the col i
    n_elem = mat->p[i+1] - mat->p[i];

    // For each element, get their column index
    for (j=0; j<n_elem; j++)
    {
      // The only difference with CSR lies here
      col = i;
      row = mat->i[n_elem_total+j];

      // Do the multiplication and cumulative sum of result on target position
      double vec_element = gsl_vector_get(vec, col);
      double mat_element = gsl_spmatrix_get(mat, row, col);
      double result_element = gsl_vector_get(result, row);

      // This works since gsl_vector_calloc initializes vector to zero, otherwise an initial set whould be done
      result_element += mat_element * vec_element;
      gsl_vector_set(result, row, result_element);
    }
    n_elem_total += n_elem;
  }
}

#elif ICC
#include "mkl.h"
void my_sparse_coo_icc(MKL_INT *rows, MKL_INT *cols, double *values, double *vec, double *result, unsigned int nnz)
{
  size_t i, tot_elem = nnz;
  int row, col;
  double mat_element, vec_element, result_element;

  for(i=0; i<tot_elem; i++)
  {
    row = rows[i];
    col = cols[i];
    mat_element = values[i];

    // Do the multiplication and cumulative sum of result on target position
    // double vec_element = gsl_vector_get(vec, col);
    vec_element = vec[col];
    // double result_element = gsl_vector_get(result, row);
    result_element = result[row];
    
    // This works since gsl_vector_calloc initializes vector to zero, otherwise an initial set whould be done
    result_element += mat_element * vec_element;
    // gsl_vector_set(result, row, result_element);
    result[row] = result_element;
  }
}
#endif
