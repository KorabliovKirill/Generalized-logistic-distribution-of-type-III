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
      MainDistribution(); // конструктор по-умолчанию
      MainDistribution(double nu, double mu, double lambda, int n); // конструктор c явным указанием параметров 
      MainDistribution(FILE* file); // конструктор из FILE
      MainDistribution(ifstream& file); // конструктор из ifstream
      MainDistribution(string path); // конструктор из пути std::string
      MainDistribution(const char* path); // конструктор из пути Си-строки

      double get_nu(); // получение параметра формы
      double get_mu(); // получение сдвига
      double get_lambda(); // получение масштаба

      void set_nu(double val); // установка парметра формы
      void set_mu(double val); // установка сдвига
      void set_lambda(double val); // установка масштаба

      double rand_value();
      double density(double x);
      double mathematicalExpectation();
      double dispersion();
      double asymmetry();
      double kurtosis();
      void plot();

      string save_params();
};