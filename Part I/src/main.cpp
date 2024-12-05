#include "distribution.hpp"
#include "mix_distribution.hpp"
#include "emp_distribution.hpp"
#include <cassert>

int distr();
int mixed();
int empir();
void testMain1();
void testMain2();
void testMix1();
void testMix2();
void testMix3();
void testEmp1();
void testEmp2();
void testEmp3();
int plot_num = 1;

int main(){
   int switcher = -1;
   bool exit = false;
   double x;

   cout << "=====================================================" << endl;
   testMain1();
   testMain2();
   testMix1();
   testMix2();
   testMix3();
   testEmp1();
   testEmp2();
   testEmp3();
   cout << "=====================================================" << endl;

   while(!exit){
      cout << "Select the type of distribution you want to work with" << endl;
      cout << "1. Main distribution" << endl;
      cout << "2. Mixed distribution" << endl;
      cout << "3. Empirical distribution" << endl;
      cout << "0. Exit" << endl;
      cin >> switcher;

      switch (switcher)
      {
      case 1:
         distr();
         break;
      case 2:
         mixed();
         break;
      case 3:
         empir();
         break;
      case 0:
         exit = true;
         break;
      default:
         cout << "Error. Try again!";
         break;
      }
   }
}

int distr(){
   int switcher;
   double x, nu, mu, lambda;

   cout << "Main distribution:" << endl;
   cout << "1. Enter params" << endl;
   cout << "0. Exit" << endl;
   cin >> switcher;

   switch (switcher)
   {
   case 1:
      cout << "Main distribution(Enter params for distribution):" << endl;

      cout << "Enter x:";
      cin >> x;

      cout << "Enter nu:";
      cin >> nu;
      if (nu < 1){
         cout << "Incorrect data (nu < 1)" << endl;
         return 0;
      } 

      cout << "Enter mu:";
      cin >> mu;

      cout << "Enter lambda:";
      cin >> lambda;
      if (lambda <= 0){
         cout << "Incorrect data (lambda <= 0)" << endl;
         return 0;
      }

      cout << "Density: " << distribution::density(nu, mu, lambda, x) << endl;
      cout << "Mathematical expectation: " << distribution::mathematicalExpectation(nu, mu, lambda) << endl;
      cout << "Dispersion: " << distribution::dispersion(nu, mu, lambda) << endl;
      cout << "Asymmetry: " << distribution::asymmetry(nu, mu, lambda) << endl;
      cout << "Kurtosis: " << distribution::kurtosis(nu, mu, lambda) << endl;
      distribution::plot(plot_num, nu, mu, lambda);
      plot_num++;

      return 1;

   case 0:
      return 0;
   }
}

int mixed(){
   int switcher;
   double x, nu1, mu1, lambda1, nu2, mu2, lambda2, p;

   cout << "Mixed distribution:" << endl;
   cout << "1. Enter params" << endl;
   cout << "0. Exit" << endl;
   cin >> switcher;

   switch (switcher)
   {
   case 1:
      cout << "Mixed distribution(Enter params for distribution):" << endl;

      cout << "Enter x:";
      cin >> x;

      cout << "Enter p:";
      cin >> p;
      if(p < 0 || p > 1){
         cout << "Incorrect data (p < 0 || p > 1)" << endl;
         return 0;
      }

      cout << "Enter nu1:";
      cin >> nu1;
      if (nu1 < 1){
         cout << "Incorrect data (nu1 < 1)" << endl;
         return 0;
      } 

      cout << "Enter mu1:";
      cin >> mu1;

      cout << "Enter lambda1:";
      cin >> lambda1;
      if (lambda1 <= 0){
         cout << "Incorrect data (lambda1 <= 0)" << endl;
         return 0;
      }

      cout << "Enter nu2:";
      cin >> nu2;
      if (nu2 < 1){
         cout << "Incorrect data (nu2 < 1)" << endl;
         return 0;
      } 

      cout << "Enter mu2:";
      cin >> mu2;

      cout << "Enter lambda2:";
      cin >> lambda2;
      if (lambda2 <= 0){
         cout << "Incorrect data (lambda2 <= 0)" << endl;
         return 0;
      }

      cout << "Density: " << mix_distribution::density(nu1, nu2, mu1, mu2, lambda1, lambda2, p, x) << endl;
      cout << "Mathematical expectation: " << mix_distribution::mathematicalExpectation(nu1, nu2, mu1, mu2, lambda1, lambda2, p) << endl;
      cout << "Dispersion: " << mix_distribution::dispersion(nu1, nu2, mu1, mu2, lambda1, lambda2, p) << endl;
      cout << "Asymmetry: " << mix_distribution::asymmetry(nu1, nu2, mu1, mu2, lambda1, lambda2, p) << endl;
      cout << "Kurtosis: " << mix_distribution::kurtosis(nu1, nu2, mu1, mu2, lambda1, lambda2, p) << endl;
      mix_distribution::plot(plot_num, nu1, nu2, mu1, mu2, lambda1, lambda2, p);
      plot_num++;

      return 1;

   case 0:
      return 0;
   }
}

int empir(){
   int switcher;
   double x, nu, mu, lambda;
   int n;
   double temp;
	vector<double> samples;

   cout << "Empirical distribution:" << endl;
   cout << "1. Enter params" << endl;
   cout << "0. Exit" << endl;
   cin >> switcher;

   switch (switcher)
   {
   case 1:
      cout << "Empirical distribution(Enter params for distribution):" << endl;

      cout << "Enter x:";
      cin >> x;

      cout << "Enter n: "; // размер выборки
      cin >> n;
      if (n < 1){
         cout << "Incorrect data (n < 1)" << endl;
         return 0;
      }

      cout << "Enter nu:";
      cin >> nu;
      if (nu < 1){
         cout << "Incorrect data (nu < 1)" << endl;
         return 0;
      } 

      cout << "Enter mu:";
      cin >> mu;

      cout << "Enter lambda:";
      cin >> lambda;
      if (lambda <= 0){
         cout << "Incorrect data (lambda <= 0)" << endl;
         return 0;
      }

      samples = emp_distribution::samples(n, nu, mu, lambda);

      temp = emp_distribution::density(x, samples);
      cout << "Density: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution::density(nu, mu, lambda, x)) << endl;

      temp = emp_distribution::mathematicalExpectation(samples);
      cout << "Mathematical expectation: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution::mathematicalExpectation(nu, mu, lambda)) << endl;

      temp = emp_distribution::dispersion(samples);
      cout << "Dispersion: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution::dispersion(nu, mu, lambda)) << endl;

      temp = emp_distribution::assymetry(samples);
      cout << "Asymmetry: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution::asymmetry(nu, mu, lambda)) << endl;

      temp = emp_distribution::kurtosis(samples);
      cout << "Kurtosis: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution::kurtosis(nu, mu, lambda)) << endl;

      emp_distribution::plot(plot_num, samples);
      plot_num++;

      return 1;

   case 0:
      return 0;
   }
}

void testMain1(){
   double x = 0, nu = 1, mu = 0, lambda = 1;
   assert(distribution::density(nu, mu, lambda, x) == 0.318);
   assert(distribution::dispersion(nu, mu, lambda) == 2.467);
   assert(distribution::kurtosis(nu, mu, lambda) == 2);
   distribution::plot(plot_num, nu, mu, lambda);
   plot_num++;
   cout << "Test 1 passed" << endl;
}

void testMain2(){
   double x = 0, nu = 1, mu = 1, lambda = 2;
   assert(distribution::density(nu, mu, lambda, x) == 0.141);
   assert(distribution::dispersion(nu, mu, lambda) == 9.869);
   assert(distribution::kurtosis(nu, mu, lambda) == 2);
   distribution::plot(plot_num, nu, mu, lambda);
   plot_num = 1;
   cout << "Test 2 passed" << endl;
}

void testMix1(){
   double x = 0, mu1 = 1, mu2 = 1, lambda1 = 2, lambda2 = 2, nu1 = 1, nu2 = 1, p = 0.7;
   assert(mix_distribution::density(nu1, nu2, mu1, mu2, lambda1, lambda2, p, x) == 0.141);
   assert(mix_distribution::dispersion(nu1, nu2, mu1, mu2, lambda1, lambda2, p) == 19.739);
   assert(mix_distribution::kurtosis(nu1, nu2, mu1, mu2, lambda1, lambda2, p) == -1.143);
   mix_distribution::plot(plot_num, nu1, nu2, mu1, mu2, lambda1, lambda2, p);
   plot_num++;
   cout << "Test 3 passed" << endl;
}

void testMix2(){
   double x = 0, mu1 = 0, mu2 = 1, lambda1 = 2, lambda2 = 1, nu1 = 1, nu2 = 1, p = 0.5;
   assert(mix_distribution::density(nu1, nu2, mu1, mu2, lambda1, lambda2, p, x) == 0.262);
   assert(mix_distribution::dispersion(nu1, nu2, mu1, mu2, lambda1, lambda2, p) == 5.184);
   assert(mix_distribution::kurtosis(nu1, nu2, mu1, mu2, lambda1, lambda2, p) == -3);
   mix_distribution::plot(plot_num, nu1, nu2, mu1, mu2, lambda1, lambda2, p);
   plot_num++;
   cout << "Test 4 passed" << endl;
}

void testMix3(){
   double x = 0, mu1 = 0, mu2 = 0, lambda1 = 1, lambda2 = 3, nu1 = 1, nu2 = 1, p = 0.5;
   assert(mix_distribution::density(nu1, nu2, mu1, mu2, lambda1, lambda2, p, x) == 0.212);
   assert(mix_distribution::dispersion(nu1, nu2, mu1, mu2, lambda1, lambda2, p) == 24.674);
   assert(mix_distribution::kurtosis(nu1, nu2, mu1, mu2, lambda1, lambda2, p) == -3);
   mix_distribution::plot(plot_num, nu1, nu2, mu1, mu2, lambda1, lambda2, p);
   plot_num = 1;
   cout << "Test 5 passed" << endl;
}

void testEmp1(){
	int ni[4] = { 100, 1000, 10000, 100000 };
	double densn[4], exp_valn[4], dispn[4], kurtn[4], asymn[4]; // массивы значений, полученных для эмпирического распределения
   double dens, exp_val, disp, kurt, asym; // теоретические значения (полученные для основного распределения)

   double x = 0, nu = 1.7, mu = 1.4, lambda = 4; // параметры основного распределения

   std::vector<double> samples;

   dens = distribution::density(x, nu, mu, lambda);
   exp_val = distribution::mathematicalExpectation(nu, mu, lambda);
   disp = distribution::dispersion(nu, mu, lambda);
   kurt = distribution::kurtosis(nu, mu, lambda);
   asym = distribution::asymmetry(nu, mu, lambda);

   for (int i = 0; i < 4; i++) {
      samples = emp_distribution::samples(ni[i], nu, mu, lambda); //генерация выборки
      // вычисление эмп характеристик
      densn[i] = emp_distribution::density(x, samples);
      exp_valn[i] = emp_distribution::mathematicalExpectation(samples);
      dispn[i] = emp_distribution::dispersion(samples);
      kurtn[i] = emp_distribution::kurtosis(samples);
      asymn[i] = emp_distribution::assymetry(samples);
      emp_distribution::plot(plot_num, samples);
      plot_num++;

      assert(densn[i] == 0.048);
      assert(abs(densn[i] - dens) == 0.053);
      assert(exp_valn[i] == 1.699);
      assert(abs(exp_valn[i] - exp_val) == 0.299);
      assert(dispn[i] == 11.161);
      assert(abs(dispn[i] - disp) == 5.604);
      assert(kurtn[i] == 3.150);
      assert(abs(kurtn[i] - kurt) == 1.767);
      assert(asymn[i] == -0.338);
      assert(abs(asymn[i] - asym) == 0.338);
	}
   cout << "Test 6 passed" << endl;
}

void testEmp2(){
   int ni[4] = { 100, 1000, 10000, 100000 };
	double densn[4], exp_valn[4], dispn[4], kurtn[4], asymn[4]; // массивы значений, полученных для эмпирического распределения
   double dens, exp_val, disp, kurt, asym; // теоретические значения (полученные для основного распределения)

   double x = 0, nu = 1.5, mu = 0, lambda = 1, p = 0.6, mu2 = 1.7, lambda2 = 3, nu2 = 2.7; // параметры основного распределения

   std::vector<double> samples;

   dens = mix_distribution::density(nu, nu2, mu, mu2, lambda, lambda2, p, x);
   exp_val = mix_distribution::mathematicalExpectation(nu, nu2, mu, mu2, lambda, lambda2, p);
   disp = mix_distribution::dispersion(nu, nu2, mu, mu2, lambda, lambda2, p);
   kurt = mix_distribution::kurtosis(nu, nu2, mu, mu2, lambda, lambda2, p);
   asym = mix_distribution::asymmetry(nu, nu2, mu, mu2, lambda, lambda2, p);

   for (int i = 0; i < 4; i++) {
      samples = emp_distribution::samples(ni[i], nu, mu, lambda); //генерация выборки
      // вычисление эмп характеристик
      densn[i] = emp_distribution::density(x, samples);
      exp_valn[i] = emp_distribution::mathematicalExpectation(samples);
      dispn[i] = emp_distribution::dispersion(samples);
      kurtn[i] = emp_distribution::kurtosis(samples);
      asymn[i] = emp_distribution::assymetry(samples);
      emp_distribution::plot(plot_num, samples);
      plot_num++;

      assert(densn[i] == 0.410);
      assert(abs(densn[i] - dens) == 0.006);
      assert(exp_valn[i] == -0.078);
      assert(abs(exp_valn[i] - exp_val) == 0.078);
      assert(dispn[i] == 1.148);
      assert(abs(dispn[i] - disp) == 0.121);
      assert(kurtn[i] == 3.552);
      assert(abs(kurtn[i] - kurt) == 2.023);
      assert(asymn[i] == -0.196);
      assert(abs(asymn[i] - asym) == 0.196);
	}
   cout << "Test 7 passed" << endl;
}

void testEmp3(){
   int n = 10000;
   int ni[4] = { 100, 1000, 10000, 100000 };
	double densn[4], exp_valn[4], dispn[4], kurtn[4], asymn[4]; // массивы значений, полученных для эмпирического распределения
   double dens, exp_val, disp, kurt, asym; // теоретические значения (полученные для основного распределения)

   double x = 0, nu = 1.5, mu = 0, lambda = 1; // параметры основного распределения

   std::vector<double> samples;

   dens = distribution::density(x, nu, mu, lambda);
   exp_val = distribution::mathematicalExpectation(nu, mu, lambda);
   disp = distribution::dispersion(nu, mu, lambda);
   kurt = distribution::kurtosis(nu, mu, lambda);
   asym = distribution::asymmetry(nu, mu, lambda);


   samples = emp_distribution::samples(n, nu, mu, lambda); //генерация выборки
   // вычисление эмп характеристик
   densn[0] = emp_distribution::density(x, samples);
   exp_valn[0] = emp_distribution::mathematicalExpectation(samples);
   dispn[0] = emp_distribution::dispersion(samples);
   kurtn[0] = emp_distribution::kurtosis(samples);
   asymn[0] = emp_distribution::assymetry(samples);
   emp_distribution::plot(plot_num, samples);
   plot_num++;

   densn[1] = emp_distribution::density(x, samples);
   exp_valn[1] = emp_distribution::mathematicalExpectation(samples);
   dispn[1] = emp_distribution::dispersion(samples);
   kurtn[1] = emp_distribution::kurtosis(samples);
   asymn[1] = emp_distribution::assymetry(samples);
   emp_distribution::plot(plot_num, samples);

   assert(densn[0] == 0.410);
   assert(abs(densn[0] - dens) == 0.006);
   assert(exp_valn[0] == -0.078);
   assert(abs(exp_valn[0] - exp_val) == 0.078);
   assert(dispn[0] == 1.148);
   assert(abs(dispn[0] - disp) == 0.121);
   assert(kurtn[0] == 3.552);
   assert(abs(kurtn[0] - kurt) == 2.023);
   assert(asymn[0] == -0.196);
   assert(abs(asymn[0] - asym) == 0.196);

   assert(densn[1] == 0.410);
   assert(abs(densn[1] - densn[0]) == 0.006);
   assert(exp_valn[1] == -0.078);
   assert(abs(exp_valn[1] - exp_valn[0]) == 0.078);
   assert(dispn[1] == 1.148);
   assert(abs(dispn[1] - dispn[0]) == 0.121);
   assert(kurtn[1] == 3.552);
   assert(abs(kurtn[1] - kurtn[0]) == 2.023);
   assert(asymn[1] == -0.196);
   assert(abs(asymn[1] - asymn[0]) == 0.196);

   cout << "Test 8 passed" << endl;
   plot_num = 100;
}