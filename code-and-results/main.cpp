// Main provides a menu to run class structures
#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "finitediffs.hpp"
#include "finitediffs2d.hpp"
#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <cmath>
#include <chrono>

using namespace std;
using namespace arma;


double I(double x);
double I_2D(double x, double y);
double I_2Dsine(double x, double y);

int main(int argc, char const *argv[]){
  Catch::Session().run();  // testing some numerical cases vs analytical results

  int method_solver;
  cout << "Press 1 to run explicit Euler \n";
  cout << "Press 2 to run implicit Euler \n";
  cout << "Press 3 to run implicit Crank Nicolson \n";
  cout << "Press 4 to run implicit Euler in 2D \n";
  cout << "Enter number:" << " ";
  cin >> method_solver;
  int dim;

  // Problem parameters 1D
  double T = 1;
  double dx = 0.01;
  double dt = dx*dx/3;
  int Lx = 1;
  double u0 = 0;
  double uN = 1;

  if (method_solver == 1){ // Forward Euler

    Explicit_Euler Solver;
    Solver.init(T,dt,Lx,dx,u0,uN);
    Solver.set_initial(I);
    vec u = Solver.solve();
  }

  if (method_solver == 2){ // Implicit Euler

    int method = 1;
    Implicit Solver;
    Solver.init(I,T,dt,Lx,dx,u0,uN,method);
    vec u = Solver.solve();
  }

  if (method_solver == 3){ // Crank-Nicolson

    int method = 2;
    Implicit Solver;
    Solver.init(I,T,dt,Lx,dx,u0,uN,method);
    vec u = Solver.solve();
  }

  if (method_solver == 4){ // ImplicitBE in 2D
    double T = 1.0;
    double h = 0.01;
    double dt = 0.1;

    int Lx = 1;
    int Ly = 1;

    // boundary conditions
    double u0x = 0;
    double uNx = 0;
    double u0y = 0;
    double uNy = 0;

    // iteration specifications
    int max_iter = 10e3;
    int threads = 1;

    Implicit_BE Solver;
    Solver.initialize(T,dt,Lx,Ly,h,u0x,uNx,u0y,uNy);
    Solver.set_initial(I_2Dsine);
    vec u = Solver.solve(max_iter,threads);
  }
  return 0;
}

double I(double x){ // zero initial conditino
  return 0.0;
}

double I_2D(double x, double y){ // gauss curve initial condition
  // assumes scaled case with x,y in [0,1]
  return 0.75*exp(-((x-0.5)/(0.2)*(x-0.5)/(0.2) + (y - 0.5)/0.2*(y - 0.5)/0.2));
}

double I_2Dsine(double x, double y){
  // saving time instead of calling M_PI
  return 0.75*sin(3.14159265359*x)*sin(3.14159265359*y);
}
