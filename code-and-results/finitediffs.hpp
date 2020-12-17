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
  vec u, u_n; // time vectors for stencil t+dt, t, t - dt
  vec t; // time vector
  double m_u0, m_uN; // boundary condition
  ofstream m_file_mesh; // cycles to file,  to get access

public:
  void initialize(double T, double dt, int Lx, double dx, double u0, double uN); // set up parameters
};

class Explicit_Euler: public Diffusion_Solver{

protected:

public:
  void init(double T, double dt, int Lx, double dx, double u0, double uN); // Set up parameters
  void set_initial(double I(double x)); // set up the inital condition
  void advance(); // for all steps
  vec solve(); // solves the system in time
  void open_mesh_to_file(ofstream&file); // open file
  void write_mesh_to_file(ofstream&file); // write solution to file

};

class Implicit: public Diffusion_Solver{ // making a class for implicit methods

protected:
  vec m_rhs; // right hand side vector in linear system - implicit methods
  vec m_a, m_b, m_c; // matrix vectors in implicit method
  double m_s;   // Fourier number
  int m_method; // choose which method to use

public:
  void init(double I(double x), double T, double dt, int Lx, double dx, double u0, double uN, int method); // if an init needed here, using initialize
  void set_initial(double I(double x));
  void set_fourier();
  void BE_setup_system(); // Need to do this for every time step!
  void CN_setup_system();
  void forward_substution();
  void backward_substition();
  void advance(); // Choose here what method to use implicit euler or CN
  vec solve();
  void open_mesh_to_file(ofstream&file); // open file
  void write_mesh_to_file(ofstream&file); // write solution to file
};

#endif
