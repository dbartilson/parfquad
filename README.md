# parfquad
Numerical quadrature of polynomials using parallel constructs in Fortran

## Features
* Equally-spaced trapezoid rule for quadrature of polynomials
* OpenMP parallelism for main part of quadrature loop
* Text (console) input of problem parameters, including:
   * Polynomial order
   * Polynomial coefficients (from 0th coefficient to nth)
   * Integration bounds (lower, upper)
   * Number of integration points (inclusive of bound points)
   * Number of OpenMP workers/threads

## Compilation
This was built using Intel Fortran compiler and OpenMP via Intel OneAPI. 

To build, please:
1. Navigate to the directory containing the `parfquad.f90` file in an Intel OneAPI command prompt for 64-bit systems. This is necessary to have all of the proper environment variables.
2. Run the following command to build the source files to an .exe:
   ```
   ifx /Qopenmp /real-size:64 /integer-size:64 parfquad.f90
   ```
4. Run the exe in a similar OneAPI command prompt. This is the easiest way for the executable to see the OpenMP libraries.
