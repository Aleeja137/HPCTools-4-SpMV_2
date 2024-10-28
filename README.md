# HPC Tools - Deliverable 3 - SpMV  
  
This code compares the performance of well established linear algebra library's functions (GSL), against our own implementation. The comparison is done both for dense and sparse matrices.  
  
The code is based on the skeleton provided by Prof. Emilio Padrón González (https://gitlab.citic.udc.es/emilio.padron/spmv). It is the third graded deliverable on the HPC Tools curse, the first one of the second block.  
  
## Compiling  
As a fresh start, a clean setup is recommended with the command:  
`make cleanall`  
  
For compilation, it is enough to execute the command:  
`make`   

## Executing  
For executing, it is enough to run the command:  
`./spmv [n] [d]`  

Where `n` is used to determine the size of a `n*n` matrix, and `d` is the density of non-zero elements for the sparse matrix.  

By default, `n=1024` and `d=0.25`.  
  
## Examples
When running the program with default settings, the following output is shown:  

```
Matriz size: 1024 x 1024 (1048576 elements)
261996 non-zero elements (24.99%)

Dense computation
----------------
Time taken by CBLAS dense computation: 1 ms
Time taken by my dense matrix-vector product: 1 ms
Result is ok!

Sparse computation
----------------
Time taken by GSL sparse computation: 0 ms
Time taken by my sparse matrix-vector product: 15 ms
Result is ok!
```