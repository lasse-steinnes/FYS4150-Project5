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

void Explicit_Euler::convergence_rate(double I(double x),int N_experiments){ // get convergence rate for FE
// make sure stability criteria is being met?
  double T = 10.0;
  double dx = 0.1;
  double dt = dx*dx/3; // let it be on stability criteria
  int Lx = 1;
  double u0 = 0;
  double uN = 1;
  int steps = 0;
  double L2;

  // initializing vectors
  vec numpoints = zeros<vec>(N_experiments); // number of time steps
  vec E = zeros<vec>(N_experiments);  // error vector
  vec h = zeros<vec>(N_experiments);  // step size dt
  vec r = zeros<vec>(N_experiments);  // convergence rate vector
  r(N_experiments-1)= 0.0;

  while (steps < N_experiments){
      L2 = 0.0;
      init(T,dt,Lx,dx,u0,uN);
      set_initial(I);
      vec u_num = solve();
      //cout << (m_x-u_num) << "\n";
      L2 = sqrt(m_dt*sum((m_x-u_num)%(m_x - u_num))); // % elementwise multiplication
      E(steps) = L2;
      h(steps) = m_dt;
      numpoints(steps) = m_Nt;
      dt = m_dt/((double) 2); // update step size
      dx =  sqrt(3*dt);
      //cout << "dx" << dx;
      steps += 1;
      }

  // get convergence rate
  for (int j = 1; j < N_experiments; j++){
      r(j) = log10(E(j)/E(j-1))/((double) log10(h(j)/h(j-1)));
    }

  cout  << "Convergence rates:\n"  << r <<"\n";
  cout << "L2-norm:\n"  << E << "\n";
  //cout << "step size:\n" << h << "\n";
}

// could add a write error to file in superclass
