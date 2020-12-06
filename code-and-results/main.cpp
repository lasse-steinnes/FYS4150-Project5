// Main provides a menu to run class structures
#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "finitediffs.hpp"
#include <iostream>
#include <armadillo>
#include <omp.h>
#include <stdio.h>

using namespace std;
using namespace arma;


double I(double x);

int main(int argc, char const *argv[]){
  Catch::Session().run();  // testing some numerical cases vs analytical results

  int method_solver;
  cout << "Press 1 to run explicit Euler \n";
  cout << "Press 2 to run implicit Euler \n";
  cout << "Press 3 to run implicit Crank Nicolson \n";
  cout << "Press 4 to calculate convergence rates \n";
  cout << "Enter number:" << " ";
  cin >> method_solver;

  if (method_solver == 1){ // Forward Euler
    Explicit_Euler Solver;
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    Solver.init(T,dt,Lx,dx,u0,uN);
    Solver.set_initial(I);
    vec u = Solver.solve();
  }

  if (method_solver == 2){ // Implicit Euler
    Implicit Solver;
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    int method = 1;
    Solver.init(T,dt,Lx,dx,u0,uN,method);
    vec u = Solver.solve();
  }

  if (method_solver == 3){ // Crank-Nicolson
    Implicit Solver;
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    int method = 2;
    Solver.init(T,dt,Lx,dx,u0,uN,method);
    vec u = Solver.solve();
  }

  if (method_solver == 4){ // Convergence rates

  }

  return 0;
}

double I(double x){
  return 0.0;
}
