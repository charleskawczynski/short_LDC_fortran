# short_LDC_fortran
This is a 200 line 3-D Lid-Driven Cavity incompressible flow solver in Fortran.
Details:

- Staggered grid with uniform spacing: dx=dy=dz=h

- The PPE is solved, to enforce div(u)=0, with red-black Gauss-Seidel method

- Time marching: Explicit Euler (1st order accurate)

- Spatial discretization: 2nd order accurate, 1st order accurate boundary treatment

- Pressure is treated purely implicitly

- One ghost cell surrounds all interior cells

- OpenMP parallelization is implemented

- Compile command: "gfortran -fopenmp -std=gnu -O3 -g main.f90 -o main"

Developed by Charlie Kawczynski. Email: charliekawczynski@gmail.com
