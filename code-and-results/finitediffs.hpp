// Setting up the superclass and child class structures
#ifndef FINITEDIFFS_HPP
#define FiNITEDIFFS_HPP

#include <armadillo>
#include <iostream>
#include <chrono>

using namespace arma;
using namespace std;

// setting up classes, public and protected variables/methods
class Diffusion_Solver{

protected:
  double m_T; // time
  double m_dt; // step size time
  double m_dx; // step size in space
  double m_dx2; // dx squared
  int m_Nt; // Number of time iterations
  int m_Nx; // number of spatial iterations
  double m_Lx; // span of dimension 1,x in [0,Lx]
  vec m_x; // spatial mesh
  vec u, u_n, u_nn; // time vectors for stencil t+dt, t, t - dt
  vec t; // time vector
  vec Ix, It; // index sets

public:
  void initialize(double T, double dt, int Lx, double dx); // set up parameters

};

class Explicit_Euler: public Diffusion_Solver{

protected:

public:
  void init(double T, double dt, int Lx, double dx); // Set up parameters
  void set_initial(double I(double x)); // set up the inital condition
  void advance(); // for all steps
  void solve(); // solves the system in time

};

class Implicit: public Diffusion_Solver{ // making a class for implicit methods

protected:

public:
  void init(); // if an init needed here, using initialize
  void forward_substution();
  void backward_substition();
  void solve_tridiag(); // use this method for the two implicit solvers
  void advance(); // Choose here what method to use implicit euler or CN
  void solve();
};

#endif
