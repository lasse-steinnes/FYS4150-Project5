#include "finitediffs.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace arma;
using namespace std;

void Implicit::init(double T, double dt, int Lx, double dx, double u0, double uN, int method){

    initialize(T, dt, Lx, dx, u0, uN); // calling the general initalizer

    m_method = method;
    double x0 = 0.0;
    m_Nx = m_Nx - 1; // redifining to number of internal grid points
    m_dxdx = m_dx*m_dx;

    m_a = zeros<vec>(m_Nx);
    m_b = zeros<vec>(m_Nx);      //allocating the matrix vectors a,b and c
    m_c = zeros<vec>(m_Nx);


    u_n =  zeros<vec>(m_Nx);        // time n, redefined to only include inner points
    u =  zeros<vec>(m_Nx);          // time n+1,redefined to only include inner points
    m_rhs = zeros<vec>(m_Nx);        //allocating the right hand side of lin. system
    m_x = zeros<vec>(m_Nx);          //allocating the grid points

    for (int i = 0; i < m_Nx; i++){
      m_x(i) = x0 + (i+1)*m_dx;       //filling the grid points from (0, 1),
    }                                //not including endpoints

    set_fourier(); // define m_s

    if (m_method == 1){
      BE_setup_system(); // implicit euler
    }
    if (m_method == 2){
      CN_setup_system();   // crank nicolson
    }
}

void Implicit::set_fourier(){ // define m_s from method used
  if (m_method == 1){
      m_s = m_dt/((double) m_dxdx); // implicit euler
  }

  if (m_method == 2){
      m_s = 0.5*m_dt/((double) m_dxdx); // crank nicolson
  }
}

void Implicit::BE_setup_system(){
    for (int i = 1; i < m_Nx-1; i++){ //filling right hand side
       m_rhs(i) = u_n(i);
    }
    // special first and last righ hand side
    // last addition is using the boundary conditions
    m_rhs(0) = u_n(0) + m_s*m_u0 ;
    m_rhs(m_Nx-1) = u_n(m_Nx-1) + m_s*m_uN;

    for (int i = 0; i < m_Nx; i++){  //filling the matrix vectors with ai, bi and ci
      m_a(i) = - m_s;      // lower diagonal vector
      m_b(i) = 1 + 2*m_s;  // diagonal vector
      m_c(i) = - m_s;      // upper diagonal vector
    }
}

void Implicit::CN_setup_system(){
    for (int i = 1; i < m_Nx-1; i++){ //filling right hand side
       m_rhs(i) = m_s*(u_n(i+1) + u_n(i-1)) + (1-2*m_s)*u_n(i);
    }
    // special first and last righ hand side
    // last addition is using the boundary conditions
    m_rhs(0) = m_s*u_n(1) + (1-2*m_s)*u_n(0) + 2*m_s*m_u0;
    m_rhs(m_Nx-1) = (1-2*m_s)*u_n(m_Nx-1) + \
                    m_s*u_n(m_Nx-2) + 2*m_s*m_uN;

    for (int i = 0; i < m_Nx; i++){  //filling the matrix vectors with ai, bi and ci
      m_a(i) = - m_s;      // lower diagonal vector
      m_b(i) = 1 + 2*m_s;  // diagonal vector
      m_c(i) = - m_s;      // upper diagonal vector
    }
}

void Implicit::forward_substution(){
  for (int i = 1; i < m_Nx; i++){
    m_b(i) = m_b(i) - (m_a(i-1)*m_c(i-1))/m_b(i-1);        //updating the main diagonal
    m_rhs(i) = m_rhs(i) - (m_rhs(i-1)*m_c(i-1))/m_b(i-1);      //updating the right hand side
  }
}

void Implicit::backward_substition(){
  u(m_Nx-1) = m_rhs(m_Nx-1)/m_b(m_Nx-1);       //giving the last element of the numerical solution
  for (int i = m_Nx-2; i >= 0; i--){
    u(i) = (m_rhs(i) - m_c(i)*u(i+1))/m_b(i);   //filling the rest of the numerical solutions
  }
}

void Implicit::advance(){
// use the forward and backward to advance in space
  forward_substution();
  backward_substition();
  u_n = u; // set u_n to be the u calculated
}

vec Implicit::solve(){
// method to advance in time and space, uses the advance method
  if (m_method == 1){       // implicit euler
    for (int n = 0; n < m_Nt;n++){
      advance();
      BE_setup_system();
    }
  }

  if (m_method == 2){     //  crank nicolson
    for (int n = 0; n < m_Nt;n++){
      advance();
      CN_setup_system();
    }
  }
  //cout << u;
  return u; // return the last updated version
}

void Implicit::convergence_rate(int N_experiments, int method){ // method to get convergence rate
// make sure stability criteria is being met
  double T = 20;
  double dx = 0.2;
  double dt = dx*dx/2;
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
      init(T,dt,Lx,dx,u0,uN,method);
      vec u_num = solve();
      L2 = sqrt(m_dt*sum((m_x-u_num)%(m_x - u_num))); // % elementwise multiplication
      E(steps) = L2;
      h(steps) = m_dt;
      numpoints(steps) = m_Nt;
      dt = m_dt/((double) 2);
      steps += 1;
      }

  // get convergence rate
  for (int j = 1; j < N_experiments; j++){
      r(j) = log10(E(j)/E(j-1))/log10(h(j)/(h(j-1)));
    }

  cout  << "Convergence rates:\n"  << r <<"\n";
  cout << "relative error:\n"  << E << "\n";
  cout << "step size:\n" << h << "\n";
}

// could also add a write error to file in superclass
void Implicit::open_mesh_to_file(ofstream&file){ // open file
  // open spin to file if true
  if (m_method == 1){
  string filename(string("./results/1D-solutions/") +  \
                  "1Dsol-Nx-" + to_string(m_Nx) + "-Nt-" + to_string(m_Nt) + "-BE.txt");
  file.open(filename);
  }

  if (m_method == 2){
    string filename(string("./results/1D-solutions/") +  \
                    "1Dsol-Nx-" + to_string(m_Nx) + "-Nt-" + to_string(m_Nt) + "-CN.txt");
    file.open(filename);
  }

  file << setprecision(8) << "T:" << m_T << " " << "dt" << m_dt \
      << " " << "dx:" << m_dx;
  file << "\n";
}

void Implicit::write_mesh_to_file(ofstream&file){ // write u to file
  file << m_u0 << "\n"; // left BC

  for (int i = 1; i < m_Nx;i++){ // inner points
    file << setprecision(15) << u(i);
    file << "\n";

  file << m_uN; // right boundary

  }
} // remember to close file
