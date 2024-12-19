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
      MainDistribution(); // конструктор по-умолчанию
      MainDistribution(double nu, double mu, double lambda, int n); // конструктор c явным указанием параметров 
      MainDistribution(FILE* file); // конструктор из FILE
      MainDistribution(ifstream& file); // конструктор из ifstream
      MainDistribution(string path); // конструктор из пути std::string
      MainDistribution(const char* path); // конструктор из пути Си-строки

      double get_nu() const; // получение параметра формы
      double get_mu() const; // получение сдвига
      double get_lambda() const; // получение масштаба
      int get_number() const override;

      void set_nu(double val); // установка парметра формы
      void set_mu(double val); // установка сдвига
      void set_lambda(double val); // установка масштаба

      double rand_value() const override;
      double density(const double& x) const override;
      double mathematicalExpectation() const override;
      double dispersion() const override;
      double asymmetry() const override;
      double kurtosis() const override;
      void plot() override; // параметры графика

      string save_params() const override;
      void load_params(const string& filename) override;
};