# FEM 2D Heat Transfer Solver (Transient & Steady-State)

## Project Overview

C++ implementation of the Finite Element Method (FEM) to solve the 2D heat conduction problem in both **steady-state** and **transient** regimes.

The model incorporates convection (Robin) and fixed temperature (Dirichlet) boundary conditions.

## Key Technical Features

* **Transient Analysis:** Solution of the heat equation using the Backward Euler / Crank-Nicolson time-stepping scheme.
* **Matrix Assembly:** Calculation and aggregation of the global stiffness matrix **[H]** (conduction + convection) and the heat capacity matrix **[C]**.
* **Numerical Integration:** Use of the **Gauss-Legendre Quadrature** method for the numerical integration of elemental matrices.
* **Geometry:** 4-node quadrilateral elements.

## Technologies
* C++ (C++17/20)
* Standard Template Library (STL)
