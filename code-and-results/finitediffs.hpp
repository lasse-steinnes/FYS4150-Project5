// Setting up the superclass and child class structures
#ifndef FINITEDIFFS_HPP
#define FiNITEDIFFS_HPP

#include <armadillo>
#include <iostream>
#include <chrono>

using namespace arma;
using namespace std;

// setting up classes, public and protected variables/methods
class Finite_Difference{

protected:

public:
  void initialize(); // set initial conditions
};

class Explicit_Euler: public Finite_Difference{

protected:

public:
  void init(); // if an init needed here, using initialize
  void advance();
  void solve();

};

class Implicit: public Finite_Difference{ // making a class for implicit methods

protected:

public:
  void init(); // if an init needed here, using initialize
  void forward_substution();
  void backward_substition();
  void solve_tridiag(); // use this method for the two implicit solvers
  void advance(); // Choose here what method to use implicit euler or CN
  void solve();
};

#endif
