#pragma once

#include "spec_func.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "interface.hpp"
#define PI 3.14

using namespace std;

class MainDistribution: public IDistribution, public IPersistent {
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

      double get_nu() const; // ��������� ��������� �����
      double get_mu() const; // ��������� ������
      double get_lambda() const; // ��������� ��������
      int get_number() const override;

      void set_nu(double val); // ��������� �������� �����
      void set_mu(double val); // ��������� ������
      void set_lambda(double val); // ��������� ��������

      double rand_value() const override;
      double density(const double& x) const override;
      double mathematicalExpectation() const override;
      double dispersion() const override;
      double asymmetry() const override;
      double kurtosis() const override;
      void plot() override; // ��������� �������

      string save_params() const override;
      void load_params(const string& filename) override;
};