gfortran -c lapblas_double_excerpts.f tools.f90 Farfield.f90
gfortran mainf.f90 lapblas_double_excerpts.o tools.o Farfield.o -o mainf
./mainf 
