// Setting up the superclass and child class structures
#ifndef FINITEDIFFS2D_HPP
#define FINITEDIFFS2D_HPP

#include <armadillo>
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace arma;
using namespace std;

class Diffusion_Solver2D{

protected:
  double m_T; // time
  double m_dt; // step size time
  double m_h; // step size in space dx = dy
  double m_hh; // dx and dy squared
  int m_Nt; // Number of time iterations
  int m_Nx; // number of spatial iterations in x
  int m_Ny; // number of iterations in y
  double m_Lx; // span of dimension 1,x in [0,Lx]
  double m_Ly; // span of dimension 2, y in [0,Ly]
  int m_k; // number of points in y dimension
  vec m_x; // spatial mesh x
  vec m_y; // spatial mesh y
  vec u, u_n, u_temp; // time vectors for stencil t+dt, t, t - dt (flattened matrix) and iterative method
  vec t; // time vector
  double m_u0x, m_uNx; // boundary condition x
  double m_u0y, m_uNy; // boundary condition y
  double m_s; // Fourier number
  ofstream m_file_mesh; // cycles to file,  to get access

  // Note: Could easily extend to having a function as boundary condition
  // but choose the easy way first


public:
  void initialize(double T, double dt, int Lx, int Ly, double h, double u0x, double uNx, double u0y, double uNy); //
};

class Implicit_BE : public Diffusion_Solver2D{

public:
  void set_initial(double I(double x, double y)); // Set the initial condition
  void jacobi_iteration_method(int max_iterations, int numthreads); // iterate and stop for given tolerance
  vec solve(int max_iterations, int numthreads);  // using jacobi iteration at each time step
  void open_mesh_to_file(ofstream&file);
  void write_mesh_to_file(ofstream&file); // write solution to file
};

#endif
