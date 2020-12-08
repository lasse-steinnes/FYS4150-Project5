// iniializer to 2D case
#include "finitediffs2d.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
using namespace arma;
using namespace std;

void Diffusion_Solver2D::initialize(double T, double dt, int Lx, int Ly, double h, double u0x, double uNx, double u0y, double uNy){
/*
Setting up necesseary parameters. Input:
- T: Time span from t = 0.
- dt: Homogenous step size in time.
- Lx: The length of the system (1D) [0,Lx]
- h : Homogenous step size in space
- u0, uN : Boundary conditions at u(0) and u(endpoint)
*/

m_T = T;
m_dt = dt;
m_Lx = Lx;
m_Ly = Ly;
m_h = h; // same spatial resolution in x and y
m_hh = h*h;
m_u0x = u0x;
m_uNx = uNx;
m_u0y = u0y;
m_uNy = uNy;
m_s = dt/m_hh;


// setting up spatial parameters //
m_Nx = round(m_Lx/((double) h)); // number of spatial iterations in x-dim
m_k = m_Nx + 1; // number of points in x dimension
m_x = linspace<vec>(0,Lx,m_Nx + 1);  // Nx + 1 points

m_Ny = round(m_Ly/((double) h)); // number of spatial iterations in y-dim
m_y = linspace<vec>(0,Ly,m_Ny + 1);

m_Nt = round(T/((double) m_dt));              // Integer number of time iterations
t = linspace<vec>(0,m_Nt*m_dt,m_Nt+1);        // time mesh

// Initialize vectors (2 dimensional mesh, as flattened array)
// x is rows and y is the columns
u_n =  zeros<vec>((m_Nx+1)*(m_Ny + 1));          // time n, i.e previous time step
u =  zeros<vec>((m_Nx+1)*(m_Ny + 1));            // time n + 1, so that t + dt
//u_temp = zeros<vec>((m_Nx +1)*(m_Ny + 1));        // temporal array for iterative method

// initial guess at u is zero. Index of u called as u(j*m_k + i) // j is along y axis (rows), and i is along the x-axis (columns)
//  i = 0,1,...,m_Nx      j = 0,1,...,m_Ny, m_k is number of points in the x direction

/* Set up boundary conditions */
for (int i = 0; i <= m_Nx; i++){ // filling in x boundaries at u0y and u0Ny
  u(0 + i) = u0y;
  u(m_Ny*m_k + i) = uNy;
  }

for (int j = 1; j < m_Ny; j++){ // filling in y boundaries at u0x and u0Nx
  u(j*m_k + 0) = u0x;
  u(j*m_k + m_Nx) = uNx;
  } // Note letting corners be defined as from filling in x boundary.
}
