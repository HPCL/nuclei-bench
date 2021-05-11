# Running

To build, see site-specific examples below. To run:
```
./eigentest <N> [header]
```
N is the dimension of the matrix, and if header = 1, it prints the header.


## Marc's initial experiments at LLNL:

```
$ module list
Currently Loaded Modules:
  1) StdEnv             (S)   3) lapack/3.9.0-xl-2020.11.12   5) openblas/0.3.13   7) magma/2.5.4   9) gcc/8.3.1
  2) hpctoolkit/2020.08       4) armadillo/10.2.0             6) cuda/11.1.1       8) essl/6.3.0   10) spectrum-mpi/rolling-release


xlc++-gpu eigentest.cpp -o eigentest -O3 -Wall -Wextra -Wpedantic -std=c++11 -I$CUDA_HOME/include/ -I$ARMAINC -I$MAGMA_INC -fopenmp -Wall -Wno-unused-function -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -L$CUDA_HOME/lib64 -lcublas -lcusolver -lcudart -L$LAPACK_DIR -llapack -lblas -L$MAGMA_LIB -lmagma_sparse -lmagma -L$ESSLLIBDIR64 -lessl
```

Results are [here](./cuda_11.3.0_465.19.01_linux.run).

## UO experimenmts

```
athena.nic.uoregon.edu:
export MAGMA_ROOT=$HOME/soft/magma-2.5.4
module load armadillo-mkl
module load cuda/11.2.0-gcc-8.3.1-v65ql5k

icpc -g   -fPIC -qopenmp -mpp -autho eigentest.cpp -o eigentest -O3 -Wall -Wextra -Wpedantic -std=c++11 \
-I$CUDA_HOME/include/ -I$ARMADILLO_DIR/include -I$MAGMA_ROOT/include -Wall -Wno-unused-function -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG  \
-L$MAGMA_ROOT/lib -Wl,-rpath,$MAGMA_ROOT/lib -lmagma_sparse -lmagma \
-L$CUDA_HOME/lib64 -Wl,-rpath,$CUDA_HOME/lib64 -L$MKLROOT/lib/intel64 -Wl,-rpath,$MKLROOT/lib/intel64 \
 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lstdc++ -lm -lcublas -lcusolver -lcudart 

gcc -g   -fPIC eigentest.cpp -o eigentest -O3 -Wall -Wextra -Wpedantic -std=c++11 -I$CUDA_HOME/include/ \
-I$ARMADILLO_DIR/include -I$MAGMA_ROOT/include -fopenmp -Wall -Wno-unused-function -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG  \
-L$MAGMA_ROOT/lib -Wl,-rpath,$MAGMA_ROOT/lib -lmagma_sparse -lmagma \
-L$CUDA_HOME/lib64 -Wl,-rpath,$CUDA_HOME/lib64 -L$MKLROOT/lib/intel64 -Wl,-rpath,$MKLROOT/lib/intel64 \
 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lstdc++ -lm -lcublas -lcusolver -lcudart 

-lmkl_intel_ilp64 instead of mkl_intel_lp64 if expecting very large arrays
```


