# Numerics101

---------------------------------------------------------------
Overview
---------------------------------------------------------------
This package is a basic C++ package which serves as a learning <br> 
tool for fundamental numerical computations techniques. It involves examples <br> 
on the use of the PETSc solver to solve block CSR matrices in serial and in parallel. <br>


---------------------------------------------------------------
Build
---------------------------------------------------------------
Provided is an examples directory for three activities. You may want to compile each using <br>
cmake. Change directory to one of the examples, say petscSolver, then execute the following commands: <br>
mkdir build <br> cd build <br> cmake .. <br> make <br>
then run the executable which will appear in the same directory <br> <br>

Note that MPI is required to compile all the examples. PETSC, however, is optional, and is <br>
only required for the petsc related examples. <br>


Dr. Mhamad Mahdi Alloush, <br> Computational Mechanics Lab, <br> American University of Beirut
