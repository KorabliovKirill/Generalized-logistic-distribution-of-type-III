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
   MainDistribution distribution;
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

      distribution.set_nu(nu);
      distribution.set_mu(mu);
      distribution.set_lambda(lambda);

      cout << "Density: " << distribution.density(x) << endl;
      cout << "Mathematical expectation: " << distribution.mathematicalExpectation() << endl;
      cout << "Dispersion: " << distribution.dispersion() << endl;
      cout << "Asymmetry: " << distribution.asymmetry() << endl;
      cout << "Kurtosis: " << distribution.kurtosis() << endl;
      distribution.plot();
      plot_num++;

      return 1;

   case 0:
      return 0;
   }
}

int mixed(){
   int switcher;
   mix_distribution::MixDistribution MD;
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

      MD.set_p(p);
      MD.get_M1().set_nu(nu1);
      MD.get_M1().set_mu(mu1);
      MD.get_M1().set_lambda(lambda1);
      MD.get_M2().set_nu(nu2);
      MD.get_M2().set_mu(mu2);
      MD.get_M2().set_lambda(lambda2);

      cout << "Density: " << MD.density(x, MD.get_M1(),MD.get_M2(), MD.get_p()) << endl;
      cout << "Mathematical expectation: " << MD.mathematicalExpectation(MD.get_M1(),MD.get_M2(), MD.get_p()) << endl;
      cout << "Dispersion: " << MD.dispersion(MD.get_M1(),MD.get_M2(), MD.get_p()) << endl;
      cout << "Asymmetry: " << MD.assymetry(MD.get_M1(),MD.get_M2(), MD.get_p()) << endl;
      cout << "Kurtosis: " << MD.kurtosis(MD.get_M1(),MD.get_M2(), MD.get_p()) << endl;
      plot_num++;

      return 1;

   case 0:
      return 0;
   }
}

int empir(){
   int switcher;
   double x, nu, mu, lambda;
   mix_distribution::MixDistribution MD;
	emp_distribution::EmpDistribution ED;
   MainDistribution distribution;
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

      distribution.set_nu(nu);
		distribution.set_mu(mu);
		distribution.set_lambda(lambda);
		ED.set_n(n);

      samples = ED.set_samples(ED.get_n(), distribution);

      temp = ED.density(x, samples);
      cout << "Density: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution.density(x)) << endl;

      temp = ED.mathematicalExpectation(samples);
      cout << "Mathematical expectation: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution.mathematicalExpectation()) << endl;

      temp = ED.dispersion(samples);
      cout << "Dispersion: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution.dispersion()) << endl;

      temp = ED.assymetry(samples);
      cout << "Asymmetry: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution.asymmetry()) << endl;

      temp = ED.kurtosis(samples);
      cout << "Kurtosis: " << temp << endl;
      cout << "Delta: " << abs(temp - distribution.kurtosis()) << endl;

      return 1;

   case 0:
      return 0;
   }
}

void testMain1(){
   MainDistribution distribution;
   double x = 0, nu = 1, mu = 0, lambda = 1;
   assert(distribution.density(x) == 0.318);
   assert(distribution.dispersion() == 2.467);
   assert(distribution.kurtosis() == 2);
   distribution.plot();
   plot_num++;
   cout << "Test 1 passed" << endl;
}

void testMain2(){
   MainDistribution distribution;
   double x = 0, nu = 1, mu = 1, lambda = 2;
   assert(distribution.density(x) == 0.141);
   assert(distribution.dispersion() == 9.869);
   assert(distribution.kurtosis() == 2);
   distribution.plot();
   plot_num = 1;
   cout << "Test 2 passed" << endl;
}

void testMix1(){
   mix_distribution::MixDistribution MD;
   double x = 0, mu1 = 1, mu2 = 1, lambda1 = 2, lambda2 = 2, nu1 = 1, nu2 = 1, p = 0.7;
   MD.set_p(p);
   MD.get_M1().set_nu(nu1);
   MD.get_M1().set_mu(mu1);
   MD.get_M1().set_lambda(lambda1);
   MD.get_M2().set_nu(nu2);
   MD.get_M2().set_mu(mu2);
   MD.get_M2().set_lambda(lambda2);
   assert(MD.density(x, MD.get_M1(), MD.get_M2(), MD.get_p()) == 0.141);
   assert(MD.dispersion(MD.get_M1(), MD.get_M2(), MD.get_p()) == 19.739);
   assert(MD.kurtosis(MD.get_M1(), MD.get_M2(), MD.get_p()) == -1.143);
   MD.plot(plot_num, MD.get_M1(), MD.get_M2(), MD.get_p());
   plot_num++;
   cout << "Test 3 passed" << endl;
}

void testMix2(){
   mix_distribution::MixDistribution MD;
   double x = 0, mu1 = 0, mu2 = 1, lambda1 = 2, lambda2 = 1, nu1 = 1, nu2 = 1, p = 0.5;
   MD.set_p(p);
   MD.get_M1().set_nu(nu1);
   MD.get_M1().set_mu(mu1);
   MD.get_M1().set_lambda(lambda1);
   MD.get_M2().set_nu(nu2);
   MD.get_M2().set_mu(mu2);
   MD.get_M2().set_lambda(lambda2);
   assert(MD.density(x, MD.get_M1(), MD.get_M2(), MD.get_p()) == 0.262);
   assert(MD.dispersion(MD.get_M1(), MD.get_M2(), MD.get_p()) == 5.184);
   assert(MD.kurtosis(MD.get_M1(), MD.get_M2(), MD.get_p()) == -3);
   MD.plot(plot_num, MD.get_M1(), MD.get_M2(), MD.get_p());
   plot_num++;
   cout << "Test 4 passed" << endl;
}

void testMix3(){
   mix_distribution::MixDistribution MD;
   double x = 0, mu1 = 0, mu2 = 0, lambda1 = 1, lambda2 = 3, nu1 = 1, nu2 = 1, p = 0.5;
   MD.set_p(p);
   MD.get_M1().set_nu(nu1);
   MD.get_M1().set_mu(mu1);
   MD.get_M1().set_lambda(lambda1);
   MD.get_M2().set_nu(nu2);
   MD.get_M2().set_mu(mu2);
   MD.get_M2().set_lambda(lambda2);
   assert(MD.density(x, MD.get_M1(), MD.get_M2(), MD.get_p()) == 0.212);
   assert(MD.dispersion(MD.get_M1(), MD.get_M2(), MD.get_p()) == 24.674);
   assert(MD.kurtosis(MD.get_M1(), MD.get_M2(), MD.get_p()) == -3);
   MD.plot(plot_num, MD.get_M1(), MD.get_M2(), MD.get_p());
   plot_num = 1;
   cout << "Test 5 passed" << endl;
}

void testEmp1(){
	int ni[4] = { 100, 1000, 10000, 100000 };
	double densn[4], exp_valn[4], dispn[4], kurtn[4], asymn[4]; // массивы значений, полученных для эмпирического распределения
   double dens, exp_val, disp, kurt, asym; // теоретические значения (полученные для основного распределения)

   mix_distribution::MixDistribution MD;
	emp_distribution::EmpDistribution ED;
   MainDistribution distribution;

   double x = 0, nu = 1.7, mu = 1.4, lambda = 4; // параметры основного распределения

   std::vector<double> samples;

   dens = distribution.density(x);
   exp_val = distribution.mathematicalExpectation();
   disp = distribution.dispersion();
   kurt = distribution.kurtosis();
   asym = distribution.asymmetry();

   for (int i = 0; i < 4; i++) {
      distribution.set_nu(nu);
		distribution.set_mu(mu);
		distribution.set_lambda(lambda);
		ED.set_n(ni[i]);
      samples = ED.set_samples(ED.get_n(), distribution); //генерация выборки
      // вычисление эмп характеристик
      densn[i] = ED.density(x, samples);
      exp_valn[i] = ED.mathematicalExpectation(samples);
      dispn[i] = ED.dispersion(samples);
      kurtn[i] = ED.kurtosis(samples);
      asymn[i] = ED.assymetry(samples);
      ED.plot(plot_num, samples);
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

   mix_distribution::MixDistribution MD;
	emp_distribution::EmpDistribution ED;
   MainDistribution distribution;

   double x = 0, nu = 1.5, mu = 0, lambda = 1, p = 0.6, mu2 = 1.7, lambda2 = 3, nu2 = 2.7; // параметры основного распределения

   std::vector<double> samples;

   MD.set_p(p);
   MD.get_M1().set_nu(nu);
   MD.get_M1().set_mu(mu);
   MD.get_M1().set_lambda(lambda);
   MD.get_M2().set_nu(nu2);
   MD.get_M2().set_mu(mu2);
   MD.get_M2().set_lambda(lambda2);

   dens = MD.density(x, MD.get_M1(),MD.get_M2(), MD.get_p());
   exp_val = MD.mathematicalExpectation(MD.get_M1(),MD.get_M2(), MD.get_p());
   disp = MD.dispersion(MD.get_M1(),MD.get_M2(), MD.get_p());
   kurt = MD.kurtosis(MD.get_M1(),MD.get_M2(), MD.get_p());
   asym = MD.assymetry(MD.get_M1(),MD.get_M2(), MD.get_p());

   for (int i = 0; i < 4; i++) {
      distribution.set_nu(nu);
		distribution.set_mu(mu);
		distribution.set_lambda(lambda);
		ED.set_n(ni[i]);
      samples = ED.set_samples(ED.get_n(), distribution); //генерация выборки
      // вычисление эмп характеристик
      densn[i] = ED.density(x, samples);
      exp_valn[i] = ED.mathematicalExpectation(samples);
      dispn[i] = ED.dispersion(samples);
      kurtn[i] = ED.kurtosis(samples);
      asymn[i] = ED.assymetry(samples);
      ED.plot(plot_num, samples);
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

   mix_distribution::MixDistribution MD;
	emp_distribution::EmpDistribution ED;
   MainDistribution distribution;

   double x = 0, nu = 1.5, mu = 0, lambda = 1; // параметры основного распределения

   std::vector<double> samples;

   dens = distribution.density(x);
   exp_val = distribution.mathematicalExpectation();
   disp = distribution.dispersion();
   kurt = distribution.kurtosis();
   asym = distribution.asymmetry();

   distribution.set_nu(nu);
   distribution.set_mu(mu);
   distribution.set_lambda(lambda);
   ED.set_n(n);
   samples = ED.set_samples(ED.get_n(), distribution); //генерация выборки
   // вычисление эмп характеристик
   densn[0] = ED.density(x, samples);
   exp_valn[0] = ED.mathematicalExpectation(samples);
   dispn[0] = ED.dispersion(samples);
   kurtn[0] = ED.kurtosis(samples);
   asymn[0] = ED.assymetry(samples);
   ED.plot(plot_num, samples);
   plot_num++;

   densn[0] = ED.density(x, samples);
   exp_valn[0] = ED.mathematicalExpectation(samples);
   dispn[0] = ED.dispersion(samples);
   kurtn[0] = ED.kurtosis(samples);
   asymn[0] = ED.assymetry(samples);
   ED.plot(plot_num, samples);
   plot_num++;

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