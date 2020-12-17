#include "catch.hpp"
#include "finitediffs.hpp"
#include "finitediffs2d.hpp"
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

TEST_CASE("Testing 1D in stationary limit"){

  /*
  test the 1D methods in the
  stationary limit u(x,t) = x.
  */

  /* Numerical results */

  auto I = [](double x) { // initial condition
        return 0.0;
      };

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
  Solver2.init(I,T,dt,Lx,dx,u0,uN,method);
  u_BE = Solver2.solve();

  // Crank-Nicolson
  method = 2;
  Implicit Solver3;
  Solver3.init(I,T,dt,Lx,dx,u0,uN,method);
  u_CN = Solver3.solve();

  /* Analytical expression in the stationary limit */
  vec u_ana = linspace<vec>(0,Lx,Nx+1);

  // Comparing analyical with numerical with tolerance
  double tol = 1E-14;
  for (int i = 0; i < Nx;i++){
    REQUIRE(abs(u_FE(i) - u_ana(i)) < tol);
  }

  for (int i = 1; i < Nx;i++){
    REQUIRE(abs(u_BE(i-1) - u_ana(i)) < tol);
    REQUIRE(abs(u_CN(i-1) - u_ana(i)) < tol);
  }
};

TEST_CASE("Testing 2D in stationary limit"){

  /*
  test the 2D Implicit (Backward) Euler method in the
  stationary limit u(x,t) = 0.
  */

  /* Numerical results */
  auto I_2D = [](double x, double y) { // initial condition
        return  0.75*exp(-((x-0.5)/(0.2)*(x-0.5)/(0.2) + (y - 0.5)/0.2*(y - 0.5)/0.2));
      };

  double T = 20;
  double h = 0.1;
  double dt = h*h/3;

  int Lx = 1;
  int Ly = 1;

  // boundary conditions
  double u0x = 0;
  double uNx = 0;
  double u0y = 0;
  double uNy = 0;

  // Get number of intervals for one dimension
  double Nh = round(Lx/((double) h));

  // iteration specifications
  int max_iter = 10e3;
  int threads = 1;

  // solve numerically
  Implicit_BE Solver;
  Solver.initialize(T,dt,Lx,Ly,h,u0x,uNx,u0y,uNy);
  Solver.set_initial(I_2D);
  vec u_num = Solver.solve(max_iter, threads); // u is a flattened matrix


  /* Analytical solution in stationary limit */
  vec u_ana = zeros<vec>((Nh+1)*(Nh+1));

  //   checking u_num = 0 in stationary limit within tolerance
  double test_tol = 1E-15;
  for (int i = 0; i < (Nh+1)*(Nh+1) ;i++){
      REQUIRE(abs(u_num(i) - u_ana(i)) < test_tol);
    }
};
