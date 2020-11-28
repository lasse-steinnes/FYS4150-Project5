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


void menu();

int main(int argc, char const *argv[]){
  menu();
  //Catch::Session().run();  // testing some numberically exact cases
  return 0;
}

void menu(){
}
