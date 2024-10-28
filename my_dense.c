#include "spmv.h"

int my_dense(const unsigned int n, const double mat[], double vec[], double result[])
{
  // code your own solver
  unsigned int i, j;
  double sum;
  for(i=0; i<n; i++)
  {
    sum = 0.0;
    for(j=0; j<n; j++)
      sum += mat[i*n+j]*vec[j];
    result[i] = sum;
  }
  return 0;
}
