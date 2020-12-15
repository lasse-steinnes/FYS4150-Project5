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

void menu();

int main(int argc, char const *argv[]){
  menu();
  //Catch::Session().run();  // testing some numberically exact cases
  int method_solver;
  cout << "Press 1 to run explicit Euler \n";
  cout << "Press 2 to run implicit Euler \n";
  cout << "Press 3 to run implicit Crank Nicolson \n";
  cout << "Enter number:" << " ";
  cin >> method_solver;
  int dim;

  if (method_solver == 1){
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    //cout << "Solve for 1 or 2 dimensions? (1/2):" << " ";
    //cin >> dim;
    /*if (dim == 1){
      cout << "er i dim1" << "\n";
      Explicit_Euler Solver;
      Solver.init(T,dt,Lx,dx,u0,uN);
      Solver.set_initial(I);
      Solver.solve();
    }*/
    cout << "er i dim2" << "\n";
    Explicit_EulerTwoDim Solver;
    Solver.init2D(T,dt,Lx,dx,u0,uN);
    Solver.set_initial2D();
    cout << "er i dim2" << "\n";
    Solver.solve2D();
  }

  if (method_solver == 2){
    Implicit Solver;
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    int method = 1;
    Solver.init(T,dt,Lx,dx,u0,uN,method);
    Solver.solve();
  }

  if (method_solver == 3){
    Implicit Solver;
    double T = 0.1;
    double dx = 0.1;
    double dt = dx*dx/3;
    int Lx = 1;
    double u0 = 0;
    double uN = 1;
    int method = 2;
    Solver.init(T,dt,Lx,dx,u0,uN,method);
    Solver.solve();
  }
  return 0;
}

void menu(){
}

double I(double x){
  return 0.0;
}
