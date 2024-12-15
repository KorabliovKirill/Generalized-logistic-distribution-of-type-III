#include "distribution.hpp"

namespace distribution{
   // метод исключения
   double rand_value(double nu, double mu, double lambda){ 
      if (nu >= 1 && nu <2){
         double r1 = 0, r2 = 0, num = 0;
         do{
            do r1 = static_cast<double>(rand()) / RAND_MAX; while (r1 == 0 || r1 == 1);
            do r2 = static_cast<double>(rand()) / RAND_MAX; while (r2 == 0 || r2 == 1);
            double x =  log(tan((PI * r1) / 2));
            num = pow(cosh(x), 1 - nu);

            if (r2 <= num)
               return mu + lambda * x;
         }while(r2 > num);
      }
      else if(nu >=2){
         double r1 = 0, r2 = 0, num = 0;
         do{
            do r1 = static_cast<double>(rand()) / RAND_MAX; while (r1 == 0 || r1 == 1);
            do r2 = static_cast<double>(rand()) / RAND_MAX; while (r2 == 0 || r2 == 1);
            double x =  log(r1/(1-r1)) / 2;
            num = pow(cosh(x), 2 - nu);

            if (r2 <= num)
               return mu + lambda * x;
         }while(r2 > num);
      }
      else{
         throw;
      }
   }

   double density(double nu, double mu, double lambda, double x){
      if (nu < 0)
         throw;
      return (1/lambda)*(1/(pow(2, nu-1)*((gamma(nu/2)*gamma(nu/2))/gamma(nu))))*(1/(pow(cosh((x-mu)/lambda), nu)));
   }

   double mathematicalExpectation(double nu, double mu, double lambda){
      if (nu < 0)
         throw;
      return mu;
   }

   double dispersion(double nu, double mu, double lambda){
      if (nu < 0)
         throw;
      return 0.5 * trigamma(nu/2) * pow(lambda, 2);
   }

   double asymmetry(double nu, double mu, double lambda){
      if (nu < 0)
         throw;
      return 0;
   }

   double kurtosis(double nu, double mu, double lambda){
      if (nu < 0)
         throw;
      return (0.5 * pentagamma(nu/2))/pow(trigamma(nu/2), 2);
   }

   void plot(int& n, double nu, double mu, double lambda){
      int N = 1000;
		double x1 = -10;
		double x2 = 10;
		double dx = (x2 - x1) / (double)N;
		std::string filename = "../../res/distribution_fplot";
		filename.append(std::to_string(n)).append(".txt");
		std::ofstream file(filename.c_str());
		for (int i = 0; i < N; i++) {

			file << x1 + i * dx << '\t' << distribution::density(nu, mu, lambda, x1 + i * dx) << std::endl;
		}
		file.close();
   }
}