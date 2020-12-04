// superclass methods
#include "finitediffs.hpp"
#include <cmath>
#include <armadillo>

void Diffusion_Solver::initialize(double T, double dt, int Lx, double dx, double u0, double uN){
  /*
  Setting up necesseary parameters. Input:
  - T: Time span from t = 0.
  - dt: Homogenous step size in time.
  - Lx: The length of the system (1D) [0,Lx]
  - dx : Homogenous step size in space
  - u0, uN : Boundary conditions at u(0) and u(endpoint)
  */

  m_T = T;
  m_dt = dt;
  m_Lx = Lx;
  m_dx = dx; //
  m_dxdx = m_dx*m_dx;
  m_u0 = u0;
  m_uN = uN;

  // setting up spatial parameters
  m_Nx = round(m_Lx/((double) m_dx)); // number of spatial iterations
  m_x = linspace<vec>(0,Lx,m_Nx+1);  // +1 because then we can use index Nx for last element

  m_Nt = round(T/((double) m_dt));              // Integer number of time iterations
  t = linspace<vec>(0,m_Nt*m_dt,m_Nt+1);              // time mesh

  // Initialize vectors (1 dimensional mesh
  //u_nn = zeros<vec>(m_Nx+1);           // time n-1, so that t -dt
  u_n =  zeros<vec>(m_Nx+1);       // time n, so that
  u =  zeros<vec>(m_Nx+1);             // time n+1, so that t + dt

  // needing to use ghost cells i = -1,0 ... Nx+1 for implicit method?

  // Index sets, only iterates over INNER points
  Ix = linspace<vec>(1,u.n_elem-1, u.n_elem-1); //(end-start +1) to get endpoint
  It = linspace<vec>(0,t.n_elem,t.n_elem + 1);  // time index set
}
