# FYS4150-Project5
Git for Project 5 in Computational Physics (FYS4150).

### Main overview
* The programs in this repository aim at solving the Diffusion equation using three numerical methods, namely Explicit Forward Euler, Implicitt Backward Euler and the Crank-Nicholson method. The description for solving the Diffuision equation is described in : [Project 5 - Solving the Diffusion equation ](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/report/Project-description-DiffusionEquation.pdf). The final report can be found at: [Aaby Rashid Steinnes - Project5](https://github.com/Seedsiz/FYS4150-Project4/blob/main/report/Aaby_Steinnes_Rashid_exploring_the_ising_model_report.pdf).

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
  5. convergence_rate: Finds the convergence rate of the explicit forward euler method.
  6. The other methods provided are write to file methods.

- [implicit.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/implicit.cpp): The subclass which solves the diffusion equation using both implicit backward euler and implicit Crank-Nicholson methods. Subclass methods provided are given in following order
  1. init: Sets up member parameters and vectors by calling the superclass method initialize from finitediffs.cpp. In addition, the method also redifines number of spatial grid points and some vectors since the implicit methods only uses the inner points of the system. An option is also made to initialize the system either for backward euler or Crank-Nicholson.
  2. BN_setup_system: Initializes the tridiagonal matrix vectors and the right hand side of the equation, when applying the backward euler method. 
  3. CN_setup_system: Initializes the tridiagonal matrix vectors and the right hand side of the equation, when applying the Crank-Nicholson method. 
  4. forward_substitution: Performing gauss elimination to the matrix equation and update the right hand side of the equation. 
  5. backward_substitution: main algorithm, finds the solutions u for a given moment in time, by performing backward substitution to the equation.
  7. advance: Moving the system to a new moment in time and calls the methods forward_substitution and backward_substitution to tolve the new solution vector u.
  8. solve: Advances the system multiple moments in time by calling the method advance Nt times. For every advancement it sets up the system by calling either BN_setup_system or CN_setup_system depending on the numerical method applied. Numerical solutions are also written to file here.
  9. The other methods provided are write to file methods.

- [finitediffs2d.hpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/finitediffs2d.hpp) Headerfile for the superclass Diffusion_Solver2D for 2 dimensions, with subclasses Implicit_BE.
- [finitediffs2d.cpp](https://github.com/lasse-steinnes/FYS4150-Project5/blob/main/code-and-results/finitediffs2d.cpp) 
- [plotspin.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plotspin.py) makes a heatmap of the final spin configuration for a given temperature.
- [plotexpvalues.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plotexpvalues.py) plots the heat capacity, magnetic susceptibility and other parameters of the system.
- [plottime.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plottime.py) makes a plot for different runtime with and without flags/parallelization.

The files can be compiled with "make all". And plots can be generated with "make argument."

How to run the programmes to reproduce the results discussed in the article: The menu provides user input on the following variables
  1. Integer number of spin particles for an axis, eg. 2,20,40,60,80,100
  2. Start point temperature, eg. 1
  3. Endpoint temperature, eg. 2.4
  4. Integer number of temperature points to be evaluated within a thread, eg. 10
  5. integer number of MC cycles, eg. 100 thousands or 1-10 million to get a good result
  6. Enter integer number of calibration cycles, 20 000 to be on the safe side.
  7. Enter integer number of threads (1 if not parallelization wanted), a value between 1-4.

After experimentation with the MC cycles, the authors decided to use 20 000 calibration cycles, to get more accurate expectation values for less number of cycles. This ensures also good results for larger spin systems and higher temperatures, which needs more calibration cycles compared to lower temperatures and smaller spin system. If one want to plot more than just expectation values, one must go in the isingmodel.cpp to decomment variables m_accepted, E_cycles and M_cycles. Then you must set the bool objects to true.

### Links and packages
- The Mersenne Twister (pseudo)random number generator was used in generating uniform distribution to draw indices and acceptance criteria. Documentation on the class mt19937_64 can be found [here.](https://www.cplusplus.com/reference/random/mt19937_64/)

- Documentation for Matplotlib from python from [here](https://matplotlib.org/)

- Documentation on seaborn displot [here](https://seaborn.pydata.org/generated/seaborn.displot.html#seaborn.displot) and heatmaps [here](https://seaborn.pydata.org/generated/seaborn.heatmap.html)

- Documentation on pandas dataframes can be found [here](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)

- Documentation on parallelization with OpenMP can be found [here](https://www.openmp.org/wp-content/uploads/OpenMP-4.5-1115-CPP-web.pdf) or for more versions [here](https://www.openmp.org/resources/refguides/)

- Documentation on armadillo can be found [here](http://arma.sourceforge.net/docs.html)
