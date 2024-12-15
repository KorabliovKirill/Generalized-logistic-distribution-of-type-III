#pragma once

#include "spec_func.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#define PI 3.14

using namespace std;

class MainDistribution {
   private:
      double nu;
      double mu;
      double lambda;
      int n;
      int numerator;

   public:
      MainDistribution(); // ����������� ��-���������
      MainDistribution(double nu, double mu, double lambda, int n); // ����������� c ����� ��������� ���������� 
      MainDistribution(FILE* file); // ����������� �� FILE
      MainDistribution(ifstream& file); // ����������� �� ifstream
      MainDistribution(string path); // ����������� �� ���� std::string
      MainDistribution(const char* path); // ����������� �� ���� ��-������

      double get_nu(); // ��������� ��������� �����
      double get_mu(); // ��������� ������
      double get_lambda(); // ��������� ��������

      void set_nu(double val); // ��������� �������� �����
      void set_mu(double val); // ��������� ������
      void set_lambda(double val); // ��������� ��������

      double rand_value();
      double density(double x);
      double mathematicalExpectation();
      double dispersion();
      double asymmetry();
      double kurtosis();
      void plot();

      string save_params();
};