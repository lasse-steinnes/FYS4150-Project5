#include "finitediffs.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace arma;
using namespace std;

void Explicit_Euler::init(double T, double dt, int Lx, double dx, double u0, double uN){
  initialize(T, dt, Lx, dx, u0, uN);
}

void Explicit_Euler::set_initial(double I(double x)){  // set up the inital condition
  for (int i = 1; i < m_Nx; i++){
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
  double dx = 0.2;
  double dt = dx*dx/3; // let it be on stability criteria
  int Lx = 1;
  double u0 = 0;
  double uN = 1;
  int steps = 0;
  double error,next_error;

  // initializing vectors
  vec numpoints = zeros<vec>(N_experiments); // number of time steps
  vec E = zeros<vec>(N_experiments);  // error vector
  vec h = zeros<vec>(N_experiments);  // step size dt
  vec r = zeros<vec>(N_experiments);  // convergence rate vector
  r(N_experiments-1)= 0.0;

  while (steps < N_experiments){
      init(T,dt,Lx,dx,u0,uN);
      set_initial(I);
      vec u_num = solve();
      //cout << (m_x-u_num) << "\n";
      // try to use max-norm instead maybe,

      error = 0.0;
      /* L2 norm
      for (int i = 0; i < m_Nx; i++){    //For all time steps
        next_error =  m_x(i) - u_num(i);   //Calulating the error.
        error = error + next_error*next_error;
        }
        */
      //Supremum norm
      for (int i = 0; i < m_Nx; i++){    //For all time steps
        next_error =  m_x(i) - u_num(i);   //Calulating the error.
        error = max(next_error,error);
      }

      E(steps) = sqrt(m_dt*error);
      h(steps) = m_dt;
      numpoints(steps) = m_Nt;
      dt = m_dt/((double) 2); // update step size
      dx =  sqrt(3*dt);//2
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

// could also add a write error to file in superclass

void Explicit_Euler::open_mesh_to_file(ofstream&file){ // open file
  // open spin to file if true
  string filename(string("./results/1D-solutions/") +  \
                  "1Dsol-Nx-" + to_string(m_Nx) + "-Nt-" + to_string(m_Nt) +  "-FE.txt");
  file.open(filename);
  file << setprecision(8) << "T:" << m_T << " " << "dt" << m_dt \
      << " " << "dx:" << m_dx;
  file << "\n";
}

void Explicit_Euler::write_mesh_to_file(ofstream&file){ // write u to file
  for (int i = 0; i <= m_Nx;i++){
    file << setprecision(15) << u(i);
    file << "\n";
  }
} // remember to close file
