// Main provides a menu to run class structures
#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "finitediffs.hpp"
#include "finitediffs2d.hpp"
#include <iostream>
#include <armadillo>
#include <omp.h>
#include <stdio.h>
#include <cmath>

using namespace std;
using namespace arma;


double I(double x);
double I_2D(double x, double y);

int main(int argc, char const *argv[]){
  //Catch::Session().run();  // testing some numerical cases vs analytical results

  int method_solver;
  cout << "Press 1 to run explicit Euler \n";
  cout << "Press 2 to run implicit Euler \n";
  cout << "Press 3 to run implicit Crank Nicolson \n";
  cout << "Press 4 to calculate convergence rates in 1D \n";
  cout << "Press 5 to run implicit Euler in 2D \n";
  cout << "Enter number:" << " ";
  cin >> method_solver;

  if (method_solver == 1){ // Forward Euler
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;

    Explicit_Euler Solver;
    Solver.init(T,dt,Lx,dx,u0,uN);
    Solver.set_initial(I);
    vec u = Solver.solve();
  }

  if (method_solver == 2){ // Implicit Euler
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    int method = 1;
    Implicit Solver;
    Solver.init(T,dt,Lx,dx,u0,uN,method);
    vec u = Solver.solve();
  }

  if (method_solver == 3){ // Crank-Nicolson
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    int method = 2;

    Implicit Solver;
    Solver.init(T,dt,Lx,dx,u0,uN,method);
    vec u = Solver.solve();
  }

  if (method_solver == 4){ // Convergence rates
    /* Could choose to have this as a test case
    convergence rates for 1D methods
    expected convergence rates FE: 1, BE: 1, CN: 2 as dt --> 0.
    */
    cout << "Calculating convergence rate for Forward Euler, \n \
            Backward Euler and Crank-Nicolson \n";

    // explicit methods
    Explicit_Euler convergence;
    //convergence.convergence_rate(I,8);

  //   implicit methods
    Implicit convergence1;
    //convergence1.convergence_rate(5,1); // BE
    convergence1.convergence_rate(5,2); // CN

  }

  if (method_solver == 5){ // do implicitBE in 2D
    double T = 0.1;
    double h = 0.1;
    double dt = h*h/3;
    // set upper limits for x and y,
    // PS: Jacobi is ensure to converge for square system
    int Lx = 1;
    int Ly = 1;

    // boundary conditions
    double u0x = 0;
    double uNx = 0;
    double u0y = 0;
    double uNy = 0;

    // iteration specifications
    int max_iter = 1e3;
    double tol = 1e-5;

    Implicit_BE Solver;
    Solver.initialize(T,dt,Lx,Ly,h,u0x,uNx,u0y,uNy);
    Solver.set_initial(I_2D);
    vec u = Solver.solve(max_iter, tol);
  }
  return 0;
}

double I(double x){
  return 0.0;
}

double I_2D(double x, double y){
  // assumes scaled case with x,y in [0,1]
  return exp(-((x-0.5)*(x-0.5) + (y - 0.5)*(y - 0.5))) ;// gauss curve
}
