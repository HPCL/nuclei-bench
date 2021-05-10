In the main program at the end of the file `test_denshf.f90`, you can set the speed of the code between 1 (slow) to 3 (fast). 
By default, it is set to 3. Also by default, the Makefile activates only OpenMP. You can also activate MPI by 
setting `USE_MPI=1` and `USE_MANYCORES=1` and setting the numbers `M_GRID` and `N_GRID` such as the Cartesian product 
`M_GRID x N_GRID = [whatever #MPI ranks you want]`. Again, this is patterned on HFODD.
