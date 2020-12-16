// Code for backward euler in 2D
#include "finitediffs2d.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>

using namespace arma;
using namespace std;

void Implicit_BE::set_initial(double I(double x, double y)){
  // Setting up the initial condition, avoiding boundaries
  for (int i =  1; i < m_Nx; i++){
    for (int j = 1; j < m_Ny; j++){
      u_n(j*m_k + i) = I(m_x(i),m_y(j)); // fill in for the zeroth state
    }
  }
 u = u_n; // needed for first  iteration
}

void Implicit_BE::jacobi_iteration_method(int max_iterations, double tol){
 // Jacobi iteration for spatial solution at a specific time
 // u_temp = zeros<vec>((m_Nx +1)*(m_Ny + 1)); // Temporal array for iterative method,
 // Setting u_temp to zero (could start random as well (might converge faster))
 double diff;
 int nn = m_Nx*m_Ny;
 int iterations = 0;
 // scheme
 while (iterations <= max_iterations){ //  && (diff > tol)
    u_temp = u; diff = 0.0;
    for (int i = 1; i < m_Nx; i++){
        for (int j = 1; j < m_Ny; j++){
            double delta_ij = u_temp((j+1)*m_k + i) + u_temp((j-1)*m_k + i) \
                              + u_temp(j*m_k + (i-1)) + u_temp(j*m_k + (i+1));
            u(j*m_k + i) = 1/((double) 1+4*m_s)*(m_s*delta_ij + u_n(j*m_k + i));

            //diff += fabs(u(j*m_k + i) - u_temp(j*m_k + i)); // accumulative difference (too small for values between 0 an 1)
            }
          }
       iterations ++;
       //diff /= ((double) nn);
     }
  //cout << diff << "\n";
}


vec Implicit_BE::solve(int max_iterations, double tol){
  // Advance in time
  for (int n = 0; n < m_Nt;n++){
      jacobi_iteration_method(max_iterations, tol);
      u_n = u; // update mesh
  }
  open_mesh_to_file(m_file_mesh);
  write_mesh_to_file(m_file_mesh);
  m_file_mesh.close();
  //cout << u;
  return u;
}


void Implicit_BE::open_mesh_to_file(ofstream&file){ // open file
  // open spin to file if true
  string filename(string("./results/2D-solutions/2Dsol-Nx-") \
                    + to_string(m_Nx) + "-Ny-" + to_string(m_Ny)  + "-Nt-" + to_string(m_Nt) + "-BE.txt");
  file.open(filename);
  file << setprecision(8) << "T:" << m_T << " " << "dt" << m_dt \
      << " " << "h-step-xy:" << m_h;
  file << "\n";
}

void Implicit_BE::write_mesh_to_file(ofstream&file){ // write u to file
  for (int j = 0; j <= m_Ny; j++){
    for (int i = 0; i <= m_Nx;i++){ // all points
    file << setprecision(15) << u(j*m_k + i) << " ";
    }
  file << "\n"; // for each row create a new line
  }
} // remember to close file
