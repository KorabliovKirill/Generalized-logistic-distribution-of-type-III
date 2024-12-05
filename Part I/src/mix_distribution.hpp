#pragma once

#include "distribution.hpp"

namespace mix_distribution{
   double rand_value(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
   double density(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p, double x);
   double mathematicalExpectation(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
   double dispersion(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
   double asymmetry(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
   double kurtosis(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
   void plot(int& n, double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
}