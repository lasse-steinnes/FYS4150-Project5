# FYS4150-Project5
Git for Project 5 in Computational Physics (FYS4150).

### Main overview
* The programs in this repository aim at solving the Diffusion equation using three numerical methods, namely Explicit Forward Euler, Implicit Backward Euler and the Crank-Nicholson method. The description for solving the Diffuision equation is described in : [Project 5 - Solving the Diffusion equation ](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/report/Project-description-DiffusionEquation.pdf). The final report can be found at: [Aaby Rashid Steinnes - Project5](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/report/Aaby-I-Rashid-S-Steinnes-L-Investigating-Finite-Difference-Schemes-of-the-Heat-Equation.pdf).

* The main challange is to solve the Diffusion equation both in 1 and 2 dimensions. The intuition of the problem can be thought of as finding the temperature gradient in a rod of length L = 1. For the 1D problem, the goal is to simulate over a significant good amount of time, so that the numerical solutions for the various methods approach a linear variation close to an analytical stationary solution of the problem. When implementing, different approaches to spatial steplength and timestep has been made, in order to study the resolution and stability of the different numerical methods.

* For 2 dimensions, the implicit backward euler method was applied. Both dimensions are equally long and the implementation of the iterative Jacobi's method was needed in order to solve the 2D problem with this particular numerical method. For every advancement in time, the program runs the Jacobi method in order to find the gradient at a particular moment in time. We used boundary conditions 0 on all boundaries in the 2 dimensional case for the implicit backward euler method.

* Textfiles and figures can be found in the folder results.

### Code: Link and description of programmes
- [main.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/main.cpp) : Runs the other programmes and provide user options through terminal.

 - [makefile](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/makefile) : Compiles and executes cpp files, and provides plot options of 1D and 2D solutions, amplifications for the various methods, 2D animation of the temperature gradient.

-  [finitediffs.hpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/finitediffs.hpp) : Headerfile for the superclass Diffusion_Solver for 1 dimension, with subclasses Explicit_Euler and Implicit.

- [finitediffs.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/finitediffs.cpp) : Provides the superclass method initialize, which initializes the system with number of spatial grid points Nx and number of time simulations Nt. It also defines the horizontal 1D vector x, and vectors for solution to the temperature gradient u(x,t).

- [explicit_euler.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/explicit_euler.cpp) : The subclass which solves the diffusion equation using the explicit forward euler method. Subclass methods provided are given in following order
  1. init: Sets up member parameters and vectors by calling the superclass method initialize from finitediffs.cpp.
  2. set_initial: Initializes the solution vector u_n with zero everywhere, except at the boundaries. The boundaries are set to u_n(0) = 0 and u_n(Nx) = 1. 
  3. advance: Moving the system to a new moment in time and calculates the new solution vector u, except at the boundaries which are fixed during the entire       simulation.
  4. solve: Advances the system multiple moments in time by calling the method advance Nt times. For every advancement the corresponding solutions are written to file.
  5. The other methods provided are write to file methods.

- [implicit.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/implicit.cpp): The subclass which solves the diffusion equation using both implicit backward euler and implicit Crank-Nicholson methods. Subclass methods provided are given in following order
  1. init: Sets up member parameters and vectors by calling the superclass method initialize from finitediffs.cpp. In addition, the method also redifines number of spatial grid points and some vectors since the implicit methods only uses the inner points of the system. An option is also made to initialize the system either for backward euler or Crank-Nicholson.
  2. BN_setup_system: Initializes the tridiagonal matrix vectors and the right hand side of the equation, when applying the backward euler method. 
  3. CN_setup_system: Initializes the tridiagonal matrix vectors and the right hand side of the equation, when applying the Crank-Nicholson method. 
  4. forward_substitution: Performing gauss elimination to the matrix equation and update the right hand side of the equation. 
  5. backward_substitution: main algorithm, finds the solutions u for a given moment in time, by performing backward substitution to the equation.
  7. advance: Moving the system to a new moment in time and calls the methods forward_substitution and backward_substitution to tolve the new solution vector u.
  8. solve: Advances the system multiple moments in time by calling the method advance Nt times. For every advancement it sets up the system by calling either BN_setup_system or CN_setup_system depending on the numerical method applied. Numerical solutions are also written to file here.
  9. The other methods provided are write to file methods.

- [finitediffs2d.hpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/finitediffs2d.hpp) Headerfile for the superclass Diffusion_Solver2D for 2 dimensions, with subclass Implicit_BE.
- [finitediffs2d.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/finitediffs2d.cpp) Provides the superclass method initialize, which initializes the system with number of spatial grid points Nx and Ny and number of time simulations Nt. It also defines the vectors x and y in the 2D plane. In addition the boundary conditions are determined for both the boundaries on x and y.
- [implicit_be_2d](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/implicit_be_2d.cpp) The subclass of the superclass Diffusion_Solver2D which solves the 2-dimensional diffusion equation using the implicit backward euler. Subclass methods provided are given in following order:
  1. set_initial: Initializes the solution vector u_n for the 2-dimensional diffusion equation with zero everywhere, expcept at the boundaries which are determined using the superclass method initialize from finitediffs2d.cpp.
  2. jacobi_iteration_method: Applying jacobi's algorithm in order to find the new solution vector u for a given moment in time. 
  3. solve: Solving the diffusion equation by calling jacobi_iteration_method Nt times and update the new solution vector u each time. 
  4. The other methods provided are write to file methods.
 - [test.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/plot_amplification.py) tests the 1D and 2D algoritmes for known stationary solutions.
 - [plot_sol.py](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/plot_sol.py) Plots the 1D solutions for all three numerical methods together with the analytical stationary state to compare. The file also plots the 2D solution for the implicit backward euler method.
- [plot_amplification.py](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/plot_amplification.py) Makes a plot of the amplification factor for all three numerical methods in 1D, along with the amplification factor for the analytical solution in 1D.
- [time_error_analysis.py](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/time_error_analysis.py) Investigates and plots the runtime and error for different spatial and temporal resolution (dx and dt) for the 1D case.


The files can be compiled with "make all". And plots can be generated with "make argument."

How to run the programmes to reproduce the results discussed in the article: The menu provides user input for which numerical method you want to apply in order to solve the diffusion equation. It also asks whether or not you want to run the program for either 1 or 2 dimensions.

### Links and packages
- Documentation for Matplotlib from python from [here](https://matplotlib.org/)

- Documentation on parallelization with OpenMP can be found [here](https://www.openmp.org/wp-content/uploads/OpenMP-4.5-1115-CPP-web.pdf) or for more versions [here](https://www.openmp.org/resources/refguides/)

- Documentation on armadillo can be found [here](http://arma.sourceforge.net/docs.html)
