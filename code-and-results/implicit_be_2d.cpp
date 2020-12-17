// Code for backward euler in 2D
#include "finitediffs2d.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <omp.h>

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

void Implicit_BE::jacobi_iteration_method(int max_iterations, int numthreads){
 // Jacobi iteration for spatial solution at a specific tim
 int nn = m_Nx*m_Ny;
 int iterations = 0;

 vector<int> i_start, i_end;
 i_start.push_back(1); i_end.push_back(round(m_Nx/((double) numthreads)));

for (int indices = 1; indices < numthreads; indices++){
  i_start.push_back(round(indices/((double) numthreads)*m_Nx));
  i_end.push_back(round((indices+1)/((double) numthreads)*m_Ny));
}

 // scheme
 while (iterations <= max_iterations){
    u_temp = u;
    int i, j, thread_q;
    double delta_ij;
    // parallelize region using chosen number of nodes
    omp_set_num_threads(numthreads);
    #pragma omp parallel default(shared) num_threads(numthreads) \
             private(thread_q)
    for (thread_q = 0; thread_q < numthreads; thread_q++){
      //printf("Thread rank: %d\n", omp_get_thread_num());
      for (i = i_start[thread_q]; i < i_end[thread_q]; i++){
        for (j = 1; j < m_Ny; j++){
            delta_ij = u_temp((j+1)*m_k + i) + u_temp((j-1)*m_k + i) \
                              + u_temp(j*m_k + (i-1)) + u_temp(j*m_k + (i+1));
            u(j*m_k + i) = 1/(1+4*m_s)*(m_s*delta_ij + u_n(j*m_k + i));
          }
        }
      } // end of parallel region
    iterations ++;
   } // end while
}


vec Implicit_BE::solve(int max_iterations, int num_threads){
  // Advance in time
  auto start = chrono::high_resolution_clock::now(); // Start timer
  for (int n = 0; n < m_Nt; n++){
      jacobi_iteration_method(max_iterations, num_threads);
      u_n = u; // update mesh
  }
  auto finish = chrono::high_resolution_clock::now(); // End timer
  open_mesh_to_file(m_file_mesh);
  write_mesh_to_file(m_file_mesh);
  chrono::duration<double, std::milli> time_ms = finish - start; //get in milliseconds
  cout << "BE (2D) took" << " " <<  time_ms.count() << " " << "milliseconds \n";
  m_file_mesh.close();
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
