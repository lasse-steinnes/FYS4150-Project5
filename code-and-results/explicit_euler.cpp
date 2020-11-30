#include "finitediffs.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>

void Explicit_Euler::init(double T, double dt, int Lx, double dx){
  initialize(T, dt, Lx, dx);
}

void Explicit_Euler::set_initial(double I(double x)){  // set up the inital condition
  // Note: x[i-Ix[0]] Is the right index.
  for (int i = 1; i < m_Nx; i++){ //Ix.back gets last element in vector
              u_n(i) = I(m_x(i));
    };

  // enforcing dirichlet boundary conditions
  u_n(0) = 0; u_n(m_Nx) = 1;
}

void Explicit_Euler::advance(){
  // advance in space
  for (int i = 1; i < m_Nx;i++){
    u(i) = u_n(i) + (m_dt/m_dx2)*(u_n(i+1) - 2*u_n(i) + u_n(i-1));
    }

    // Enforce Dirichlet boundary condition
    u(0) = 0; u(m_Nx) =  1;
    u_n = u; // Set u_m = u to be used in next time step (copy u over to u_n)
}

void Explicit_Euler::solve(){ // solves the system in time
  // Uses advance to solve the system for inner mesh points in time
  for (int n = 0; n < m_Nt;n++){
      advance();
  }
}
