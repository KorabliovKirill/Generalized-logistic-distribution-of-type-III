#pragma once

#include "spec_func.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#define PI 3.14

using namespace std;

namespace distribution{
   double rand_value(double nu, double mu, double lambda);
   double density(double nu, double mu, double lambda, double x);
   double mathematicalExpectation(double nu, double mu, double lambda);
   double dispersion(double nu, double mu, double lambda);
   double asymmetry(double nu, double mu, double lambda);
   double kurtosis(double nu, double mu, double lambda);
   void plot(int& n, double nu, double mu, double lambda);
}