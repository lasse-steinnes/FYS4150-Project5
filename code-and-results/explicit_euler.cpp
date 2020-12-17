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
  open_mesh_to_file(m_file_mesh);
  auto start = chrono::high_resolution_clock::now(); // Start timer

  for (int n = 0; n < m_Nt;n++){
      advance();
      write_mesh_to_file(m_file_mesh);
  }

  auto finish = chrono::high_resolution_clock::now(); // End timer
  chrono::duration<double, std::milli> time_ms = finish - start; //get in milliseconds
  cout << "FE took"  << " " << time_ms.count() << " " << "milliseconds \n";
  m_file_mesh.close();
  return u;
}


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
  for (int i = 0; i <= m_Nx; i++){
      file << setprecision(15) << u(i) << " ";
  }
  file << "\n";

} // remember to close file
