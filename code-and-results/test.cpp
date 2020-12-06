#include "catch.hpp"
#include "finitediffs.hpp"
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

TEST_CASE("Testing 1D in stationary limit") {

  auto I = [](double x) {
        return 0.0;
      };
  /*
  test the 1D methods in the
  stationary limit u(x,t) = x.
  */

  /* Numerical results */

  double T = 20;
  double dx = 0.1;
  double dt = dx*dx/3;
  int Lx = 1;
  double u0 = 0;
  double uN = 1;
  double Nx = round(Lx/((double) dx)); // number of spatial iterations
  vec u_FE, u_BE, u_CN; // numerical results
  int method;

  // Forward euler
  Explicit_Euler Solver;
  Solver.init(T,dt,Lx,dx,u0,uN);
  Solver.set_initial(I);
  u_FE = Solver.solve();

  // Backward euler
  method = 1;
  Implicit Solver2;
  Solver2.init(T,dt,Lx,dx,u0,uN,method);
  u_BE = Solver2.solve();

  // Crank-Nicolson
  method = 2;
  Implicit Solver3;
  Solver3.init(T,dt,Lx,dx,u0,uN,method);
  u_CN = Solver3.solve();

  /* Analytical expressions */

  vec u_ana = linspace<vec>(0,Lx,Nx+1);  // analytical mesh point

  // Comparing analyical with numerical with tolerance
  // Is there a better way to do this? Catch doesnt provide array
  // testing with tolerance (?)
  double tol = 1E-15;
  for (int i = 0; i < Nx;i++){
    REQUIRE(abs(u_FE(i) - u_ana(i)) < tol);
  };
  for (int i = 1; i < Nx;i++){
    REQUIRE(abs(u_BE(i-1) - u_ana(i)) < tol);
    REQUIRE(abs(u_CN(i-1) - u_ana(i)) < tol);
  }
}

TEST_CASE("Testing convergence rates 1D"){
  /*
  Test convergence rates for 1D methods
  Expected convergence rates FE: 1, BE: 1, CN: 2 as dt --> 0.
  */

}
