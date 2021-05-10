```
$ module list
Currently Loaded Modules:
  1) StdEnv             (S)   3) lapack/3.9.0-xl-2020.11.12   5) openblas/0.3.13   7) magma/2.5.4   9) gcc/8.3.1
  2) hpctoolkit/2020.08       4) armadillo/10.2.0             6) cuda/11.1.1       8) essl/6.3.0   10) spectrum-mpi/rolling-release

xlc++-gpu eigentest.cpp -o eigentest -O3 -Wall -Wextra -Wpedantic -std=c++11 -I$CUDA_HOME/include/ -I$ARMAINC -I$MAGMA_INC -fopenmp -Wall -Wno-unused-function -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -L$CUDA_HOME/lib64 -lcublas -lcusolver -lcudart -L$LAPACK_DIR -llapack -lblas -L$MAGMA_LIB -lmagma_sparse -lmagma -L$ESSLLIBDIR64 -lessl
```
