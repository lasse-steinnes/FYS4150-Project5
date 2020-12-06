#include "finitediffs.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>

using namespace arma;
using namespace std;

void Explicit_Euler::init(double T, double dt, int Lx, double dx, double u0, double uN){
  initialize(T, dt, Lx, dx, u0, uN);
}

void Explicit_Euler::set_initial(double I(double x)){  // set up the inital condition
  // Note: x[i-Ix[0]] Is the right index.
  for (int i = 1; i < m_Nx; i++){ //Ix.back gets last element in vector
              u_n(i) = I(m_x(i));
    };

  // enforcing dirichlet boundary conditions
  u_n(0) = m_u0; u_n(m_Nx) = m_uN;
}

void Explicit_Euler::advance(){
  // advance in space
  for (int i = 1; i < m_Nx;i++){
    u(i) = u_n(i) + (m_dt/m_dxdx)*(u_n(i+1) - 2*u_n(i) + u_n(i-1));
    }

    // Enforce Dirichlet boundary condition
    u(0) = m_u0; u(m_Nx) =  m_uN;
    u_n = u; // Set u_m = u to be used in next time step (copy u over to u_n)
}

vec Explicit_Euler::solve(){ // solves the system in time
  // Uses advance to solve the system for inner mesh points in time
  for (int n = 0; n < m_Nt;n++){
      advance();
  }
  //cout << u_n;
  return u;
}
