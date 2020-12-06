#include "finitediffs.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
using namespace arma;
using namespace std;

void Implicit::init(double T, double dt, int Lx, double dx, double u0, double uN, int method){

    initialize(T, dt, Lx, dx, u0, uN); // calling the general initalizer


    double x0 = 0.0;
    m_Nx = m_Nx; // redifining to number of internal grid points

    m_a = zeros<vec>(m_Nx+1);
    m_b = zeros<vec>(m_Nx+1);      //allocating the matrix vectors a,b and c
    m_c = zeros<vec>(m_Nx+1);


    /*u_n =  zeros<vec>(m_Nx+1);        // time n, redefined to only include inner points
    u =  zeros<vec>(m_Nx+1);          // time n+1,redefined to only include inner points
    m_rhs = zeros<vec>(m_Nx+1);        //allocating the right hand side of lin. system
    m_x = zeros<vec>(m_Nx+1);*/          //allocating the grid points

    /*for (int i = 0; i < m_Nx+1; i++){
      m_x(i) = x0 + (i+1)*m_dx;       //filling the grid points from 0 to 1
    }*/
    m_rhs = zeros<vec>(m_Nx+1);
    u_n(0) = 0; u_n(m_Nx) = 1;
    u(0) = 0; u(m_Nx) = 1;


    if (method == 1){ //implicit euler
      m_s = m_dt/((double) m_dxdx);

      for (int i = 1; i < m_Nx-1; i++){ //filling right hand side
         m_rhs(i) = u_n(i);
      }
      // special first and last righ hand side
      // last addition is using the boundary conditions
      m_rhs(1) = u_n(0) + m_s*m_u0 ;
      m_rhs(m_Nx-1) = u_n(m_Nx-1) + m_s*m_uN;

    }

    if (method == 2){ // crank nicolson
      m_s = 0.5*m_dt/((double) m_dxdx);

      for (int i = 1; i < m_Nx-1; i++){ //filling right hand side
         m_rhs(i) = m_s*(u_n(i+1) + u_n(i-1)) + (1-2*m_s)*u_n(i);
      }
      // special first and last righ hand side
      // last addition is using the boundary conditions
      m_rhs(0) = m_s*u_n(1) - (1-2*m_s)*u_n(0) + 2*m_s*m_u0;
      m_rhs(m_Nx-1) = (1-2*m_s)*u_n(m_Nx-1) + \
                      m_s*u_n(m_Nx-2) + 2*m_s*m_uN;
    }

    for (int i = 0; i < m_Nx; i++){  //filling the matrix vectors with ai, bi and ci
      m_a(i) = - m_s;      // lower diagonal vector
      m_b(i) = 1 + 2*m_s;  // diagonal vector
      m_c(i) = - m_s;      // upper diagonal vector
    }
}

void Implicit::forward_substution(){
  for (int i = 1; i < m_Nx; i++){
    m_b(i) = m_b(i) - (m_a(i-1)*m_c(i-1))/m_b(i-1);        //updating the main diagonal
    m_rhs(i) = m_rhs(i) - (m_rhs(i-1)*m_c(i-1))/m_b(i-1);      //updating the right hand side g_i
  }
}

void Implicit::backward_substition(){
  u(m_Nx-1) = m_rhs(m_Nx-1)/m_b(m_Nx-1);       //giving the last element of the numerical solution
  for (int i = m_Nx-2; i > 0; i--){
    u(i) = (m_rhs(i) - m_c(i)*u_n(i+1))/m_b(i);   //filling the rest of the numerical solutions
  }
  u_n = u;
}

void Implicit::advance(){
// use the forward and backward to advance in space
  forward_substution();
  backward_substition();
}

void Implicit::solve(){
// method to advance in time and space, uses the advance method
for (int n = 0; n < m_Nt;n++){
    advance();
    cout << u << "\n";
  }
}

/* Using General algo, so that we can use input which is different for the two
implicit methods (not equal rhs)


//Initializing the system for the general algorithm
void Project1::gen_Initialize(double x0, double xn, int N, double f(double x)){
  m_N = N;                              //Number of grid points
  m_h = ((double)(xn - x0)/(m_N+1));    // steplength h

  m_a = new double[m_N];
  m_b = new double[m_N];      //allocating the matrix vectors a,b and c
  m_c = new double[m_N];

  m_v = new double[m_N];     //allocating the numerical solutions v
  m_g = new double[m_N];     //allocating the right hand side of the equation g_i = h^2f_i
  m_x = new double[m_N];     //allocating the grid points
  double hh = m_h*m_h;
  for (int i = 0; i < m_N; i++){
    m_x[i] = x0 + (i+1)*m_h;       //filling the grid points from 0 to 1
    m_g[i] = hh*f(m_x[i]);         //filling g_i, using the function f_i for different x_i
  }
}

//setting the matrix elements for the general algorithm
void Project1::set_matrix_elements(double ai, double bi, double ci){
  for (int i = 0; i < m_N; i++){
    m_a[i] = ai;
    m_b[i] = bi;         //filling the matrix vectors with ai, bi and ci
    m_c[i] = ci;
  }


//setting the matrix elements for the general algorithm
void Project1::set_matrix_elements(double ai, double bi, double ci){
  for (int i = 0; i < m_N; i++){
    m_a[i] = ai;
    m_b[i] = bi;         //filling the matrix vectors with ai, bi and ci
    m_c[i] = ci;
  }
}

//general forward substitution
void Project1::gen_forward_sub(){
  for (int i = 1; i < m_N; i++){
    m_b[i] = m_b[i] - (m_a[i-1]*m_c[i-1])/m_b[i-1];      //updating the main diagonal
    m_g[i] = m_g[i] - (m_g[i-1]*m_c[i-1])/m_b[i-1];      //updating the right hand side g_i
  }
}

//general backward substitution
void Project1::gen_backward_sub(){
  m_v[m_N-1] = m_g[m_N-1]/m_b[m_N-1];       //giving the last element of the numerical solution
  for (int i = m_N-2; i >= 0; i--){
    m_v[i] = (m_g[i] - m_c[i]*m_v[i+1])/m_b[i];   //filling the rest of the numerical solutions
  }
}
*/
