#include "mix_distribution.hpp"

namespace mix_distribution{
   double rand_value(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      if(nu1 < 0 || nu2 < 0 || p < 0 || p > 1)
         throw;
      double r;
      do r = (double)rand() / RAND_MAX; while (r == 0. || r == 1.);

      if(p < r) return distribution::rand_value(nu1, mu1, lambda1);
      else return distribution::rand_value(nu2, mu2, lambda2);
   }

   double density(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p, double x){
      if(nu1 < 0 || nu2 < 0 || p < 0 || p > 1)
         throw;
      return (1 - p) * distribution::density(nu1, mu1, lambda1, x) + p * distribution::density(nu2, mu2, lambda2, x);
   }

   double mathematicalExpectation(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      if(nu1 < 0 || nu2 < 0 || p < 0 || p > 1)
         throw;
      return (1 - p) * distribution::mathematicalExpectation(nu1, mu1, lambda1) + p * distribution::mathematicalExpectation(nu2, mu2, lambda2);
   }

   double dispersion(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      if(nu1 < 0 || nu2 < 0 || p < 0 || p > 1)
         throw;

      double D1 = (1 - p) * pow(distribution::mathematicalExpectation(nu1, mu1, lambda1), 2) + distribution::dispersion(nu1, mu1, lambda1);
      double D2 = p * pow(distribution::mathematicalExpectation(nu2, mu2, lambda2), 2) + distribution::dispersion(nu2, mu2, lambda2);

      return D1 + D2 - pow(mathematicalExpectation(nu1, nu2, mu1, mu2, lambda1, lambda2, p), 2);
   }

   double asymmetry(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      if(nu1 < 0 || nu2 < 0 || p < 0 || p > 1)
         throw;

      double M1 = distribution::mathematicalExpectation(nu1, mu1, lambda1);
      double M2 = distribution::mathematicalExpectation(nu2, mu2, lambda2);
      double M = mathematicalExpectation(nu1, nu2, mu1, mu2, lambda1, lambda2, p);
      double D1 = distribution::dispersion(nu1, mu1, lambda1);
      double D2 = distribution::dispersion(nu2, mu2, lambda2);
      double el1 = (1 - p) * (pow(M1 - M, 3)) + 3 * (M1 - M) * D1 + pow(D1, 3.0 / 2.0) * distribution::asymmetry(nu1, mu1, lambda1);
      double el2 = p * (pow(M2 - M, 3)) + 3 * (M2 - M) * D2 + pow(D2, 3.0 / 2.0) * distribution::asymmetry(nu2, mu2, lambda2);

      return (el1 + el2) / pow(dispersion(nu1, mu1, lambda1, nu2, mu2, lambda2, p), 2);
   }

   double kurtosis(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      if(nu1 < 0 || nu2 < 0 || p < 0 || p > 1)
         throw;

      double M1 = distribution::mathematicalExpectation(nu1, mu1, lambda1);
      double M2 = distribution::mathematicalExpectation(nu2, mu2, lambda2);
      double M = mathematicalExpectation(nu1, nu2, mu1, mu2, lambda1, lambda2, p);
      double D1 = distribution::dispersion(nu1, mu1, lambda1);
      double D2 = distribution::dispersion(nu2, mu2, lambda2);
      double el1 = (1 - p) * (pow(M1 - M, 4) + 6 * pow(M1 - M, 2) * D1 + 4 * (M1 - M) * pow(D1, 3.0 / 2.0) * distribution::asymmetry(nu1, mu1, lambda1) + pow(D1, 2) * (distribution::asymmetry(nu1, mu1, lambda1) + 3));
		double el2 = p * (pow(M2 - M, 4) + 6 * pow(M2 - M, 2) * D2 + 4 * (M2 - M) * pow(D2, 3.0 / 2.0) * distribution::asymmetry(nu2, mu2, lambda2) + pow(D2, 2) * (distribution::asymmetry(nu2, mu2, lambda2) + 3));

		return (el1 + el2) / pow(dispersion(nu1, mu1, lambda1, nu2, mu2, lambda2, p), 2) - 3;

   }

   void plot(int& n, double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      int N = 1000;
		double x1 = -10;
		double x2 = 10;
		double dx = (x2 - x1) / (double)N;
		std::string filename = "D:/WorkSpace/C++/Projects/Generalized-logistic-distribution-of-type-III/res/mix_distribution_fplot";
		filename.append(std::to_string(n)).append(".txt");
		std::ofstream file(filename.c_str());

      if (!file.is_open()) {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return;
      }

		for (int i = 0; i < N; i++) {

			file << x1 + i * dx << '\t' << mix_distribution::density(nu1, nu2, mu1, mu2, lambda1, lambda2, p, x1 + i * dx) << std::endl;
		}
		file.close();
   }
}