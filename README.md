# HPC Tools - Deliverable 4 - SpMV  
  
This code compares the performance of well established linear algebra library's functions (GSL and MKL) against our own implementation. The comparison is done both for dense and sparse matrices.  
  
The code is based on the skeleton provided by Prof. Emilio Padrón González (https://gitlab.citic.udc.es/emilio.padron/spmv). It is the fourth graded deliverable on the HPC Tools curse, the second one of the second block.  
  
## Compiling  
As a fresh start, a clean setup is recommended with the command:  
`make cleanall`  
  
For compilation, it is enough to execute the command (by default GCC is used with no optimizations):  
`make`   

For specific compiler optimization options, and/or using different compilers, it is enough to execute one of the following options:  
```
make gcc-O0  
make gcc-O2-novec  
make gcc-O3-vec  
make gcc-Ofast-vec  
make icc-O0  
make icc-O2-novec  
make icc-O3-vec  
make icc-Ofast-vec (not working)
```  


## Executing  
For executing, it is enough to run the command:  
`./spmv [n] [d]`  

Where `n` is used to determine the size of a `n*n` matrix, and `d` is the density of non-zero elements for the sparse matrix.  

For this test, `n=16384` and `d=0.10`.  

## Results
The following time measurements were taken as the average of 5 distinct executions.  

| GCC (ms) | O0   | O2-novec | O3-vec | Ofast-vec | Ref |
| -------- | ---- | -------- | ------ | --------- | --- |
| my_dense | 735 | 344      | 343    | 342       | 342 |
| my_coo   | 4707 | 4515     | 4508   | 4485      | 74  |
| my_csr   | 7337 | 7190     | 7200   | 7170      | 35  |
| my_csc   | 7361 | 7219     | 7223   | 7216      | 32  |

| ICC (ms) | O0  | O2-novec | O3-vec | fast | Ref |
| -------- | --- | -------- | ------ | ---- | --- |
| my_dense | 737 | 173      | 160    | x    | 227 |
| my_coo   | 168 | 74       | 73     | x    | 42  |
| my_csr   | 105 | 34       | 34     | x    | 23  |
| my_csc   | 121 | 73       | 73     | x    | 73  |

It appears that the gcc compiler, aside from the my_dense function, does not optimize well my code. Maybe the other functions are not simple enough (could be rewritten in a more straightforward way) or maybe it requires help with flags. Nevertheless, the sheer sparse computation speed using GSL library impressed me.  

On the other hand, the icc compiler seems to optimize really well both MKL's and my implementation. This makes sense since the processor on FT3 computing nodes is from Intel (Intel Xeon Platinum 8352Y), and their compiler is designed to make the most of Intel's architecture knowledge about their processors. Also, for this case, it looks like my implementation was optimized a lot. The sparse computation speed is even more impressing.    

I could not make the icc-fast target work. It is strange since I used almost the same flags as the other targets, but an error pops up with `unresolved mkl_blas_lp64_daxpy`. I tried linking libraries, loading other modules, etc. Still did not work.    

## Examples
When running the program with default settings, the following output is shown for when using the gcc compiler:  

```
Matriz size: 16384 x 16384 (268435456 elements)
26837519 non-zero elements (10.00%)
Compiled with gcc compiler

Dense computation GSL
----------------
Time taken by CBLAS dense computation: 342 ms
Time taken by my dense matrix-vector product: 735 ms
Result is ok!
Sparse computation GSL
----------------

Sparse computation
----------------
Time taken by GSL COO sparse computation: 72 ms
Time taken by GSL CSR sparse computation: 34 ms
Time taken by GSL CSC sparse computation: 38 ms
COO result is ok! - Time taken by my COO sparse matrix-vector product: 4667 ms
CSR result is ok! - Time taken by my CSR sparse matrix-vector product: 7344 ms
CSC result is ok! - Time taken by my CSC sparse matrix-vector product: 7391 ms
```  

And the following output for when using the icc compiler:  
```
Matriz size: 16384 x 16384 (268435456 elements)
26837519 non-zero elements (10.00%)
Compiled with icc compiler

Dense computation MKL
----------------
Time taken by MLK BLAS dense computation: 274 ms
Time taken by my dense matrix-vector product: 161 ms
Result is ok!
Sparse computation MKL
----------------
Time taken by MKL COO sparse computation: 42 ms
Time taken by MKL CSR sparse computation: 23 ms
Time taken by MKL CSC sparse computation: 73 ms
COO result is ok! - Time taken by my COO sparse matrix-vector product: 72 ms
CSR result is ok! - Time taken by my CSR sparse matrix-vector product: 35 ms
CSC result is ok! - Time taken by my CSC sparse matrix-vector product: 82 ms
```