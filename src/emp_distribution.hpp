#pragma once

#include "distribution.hpp"
#include "mix_distribution.hpp"
#include <algorithm>

namespace emp_distribution{
   vector<double> samples(int n, double nu, double mu, double lambda);
	vector<double> samples(int n, double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p);
	vector<double> samples(int n, vector<double>& samples);
	vector<double> rand_values(int n, vector<double>& samples);
	double density(double x, vector<double>& samples);
	double mathematicalExpectation(vector<double>& samples);
	double dispersion(vector<double>& samples);
	double assymetry(vector<double>& samples);
	double kurtosis(vector<double>& samples);
	void plot(int& n, vector<double>& samples);
}