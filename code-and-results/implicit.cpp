#include "finitediffs.hpp"
#include <armadillo>
#include <iostream>
#include <chrono>
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

    setup_system();
}

void Implicit::setup_system(){
  if (m_method == 1){ //implicit euler
    m_s = m_dt/((double) m_dxdx);

    for (int i = 1; i < m_Nx-1; i++){ //filling right hand side
       m_rhs(i) = u_n(i);
    }
    // special first and last righ hand side
    // last addition is using the boundary conditions
    m_rhs(0) = u_n(0) + m_s*m_u0 ;
    m_rhs(m_Nx-1) = u_n(m_Nx-1) + m_s*m_uN;

  }

  if (m_method == 2){ // crank nicolson
    m_s = 0.5*m_dt/((double) m_dxdx);

    for (int i = 1; i < m_Nx-1; i++){ //filling right hand side
       m_rhs(i) = m_s*(u_n(i+1) + u_n(i-1)) + (1-2*m_s)*u_n(i);
    }
    // special first and last righ hand side
    // last addition is using the boundary conditions
    m_rhs(0) = m_s*u_n(1) + (1-2*m_s)*u_n(0) + 2*m_s*m_u0;
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

void Implicit::solve(){
// method to advance in time and space, uses the advance method
for (int n = 0; n < m_Nt;n++){
    advance();
    setup_system();
  }
  //cout << u;
}
