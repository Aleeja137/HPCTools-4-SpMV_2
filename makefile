CC=gcc
CFLAGS=-O0 -DGCC -Wall -Wextra
LDLIBS=-lgsl -lgslcblas

SRC=spmv.c my_dense.c my_sparse.c timer.c
OBJ=$(SRC:.c=.o)

spmv: $(OBJ)
	@echo "Compiler: $(CC)"
	@echo "Flags: $(CFLAGS)"
	$(CC) $(CFLAGS) $(LDLIBS) -o $@ $^

gcc-O0: CFLAGS = -DGCC -O0 -Wall -Wextra
gcc-O0: spmv

gcc-O2-novec: CFLAGS = -DGCC -O2 -Wall -Wextra -fno-tree-vectorize
gcc-O2-novec: spmv

gcc-O3-vec: CFLAGS = -DGCC -O3 -Wall -Wextra -ftree-loop-vectorize -ftree-slp-vectorize -march=native
gcc-O3-vec: spmv

gcc-Ofast-vec: CFLAGS = -DGCC -Wall -Wextra -Ofast -ftree-loop-vectorize -ftree-slp-vectorize -march=native
gcc-Ofast-vec: spmv

icc-O0: CFLAGS = -DICC -O0 -Wall -Wextra -diag-disable=10441
icc-O0: LDLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
icc-O0: spmv

icc-O2-novec: CFLAGS = -DICC -O2 -Wall -Wextra -qno-vectorize -diag-disable=10441
icc-O2-novec: LDLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
icc-O2-novec: spmv

icc-O3-vec: CFLAGS = -DICC -O3 -Wall -Wextra -xHost -diag-disable=10441
icc-O3-vec: LDLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
icc-O3-vec: spmv

icc-Ofast-vec: CFLAGS = -DICC -O3 -fast -Wall -Wextra  -xHost -diag-disable=10441
icc-Ofast-vec: LDLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm
icc-Ofast-vec: spmv

clean:
	$(RM) $(OBJ) *~

cleanall:
	$(RM) $(OBJ) spmv *~
