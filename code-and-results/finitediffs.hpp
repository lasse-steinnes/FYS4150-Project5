// Setting up the superclass and child class structures
#ifndef FINITEDIFFS_HPP
#define FINITEDIFFS_HPP

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
  double m_dxdx; // dx squared
  int m_Nt; // Number of time iterations
  int m_Nx; // number of spatial iterations
  double m_Lx; // span of dimension 1,x in [0,Lx]
  vec m_x; // spatial mesh
  vec u, u_n, u_nn; // time vectors for stencil t+dt, t, t - dt
  vec t; // time vector
  vec Ix, It; // index sets
  double m_u0, m_uN; // boundary condition

public:
  void initialize(double T, double dt, int Lx, double dx, double u0, double uN); // set up parameters

};

class Explicit_Euler: public Diffusion_Solver{

protected:

public:
  void init(double T, double dt, int Lx, double dx, double u0, double uN); // Set up parameters
  void set_initial(double I(double x)); // set up the inital condition
  void advance(); // for all steps
  void solve(); // solves the system in time

};

class Implicit: public Diffusion_Solver{ // making a class for implicit methods

protected:
  vec m_rhs; // right hand side vector in linear system - implicit methods
  vec m_a, m_b, m_c; // matrix vectors in implicit method
  double m_s;   // Fourier number
  int m_method; // choose which method to use

public:
  void init(double T, double dt, int Lx, double dx, double u0, double uN, int method); // if an init needed here, using initialize
  void setup_system(); // Need to do this for every time step!
  void forward_substution();
  void backward_substition();
  void advance(); // Choose here what method to use implicit euler or CN
  void solve();
};

#endif
