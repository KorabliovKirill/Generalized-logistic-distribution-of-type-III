#include "distribution.hpp"

MainDistribution::MainDistribution() : // конструктор по-умолчанию
	nu(1), mu(0), lambda(1), n(0), numerator(0)
{}

MainDistribution::MainDistribution(double nu, double mu, double lambda, int n) : // конструктор c явным указанием параметров 
	nu(nu), mu(mu), lambda(lambda), n(n), numerator(0)
{}

MainDistribution::MainDistribution(FILE* file) { // конструктор из FILE
	if (file != NULL) {
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
	}
	numerator = 0;
}

MainDistribution::MainDistribution(string path) { // конструктор из пути std:: string
	FILE* file;
	fopen_s(&file, path.c_str(), "r");
	if (file != NULL) {
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
	}
	numerator = 0;
}

MainDistribution::MainDistribution(const char* path) { // конструктор из пути Си-строки
	FILE* file;
	fopen_s(&file, path, "r");
	if (file != NULL) {
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
	}
	numerator = 0;
}

MainDistribution::MainDistribution(std::ifstream& file) { // конструктор из ifstream
	file >> nu >> mu >> lambda >> n;
	numerator = 0;
}

double MainDistribution::get_nu() const { // получение параметра формы
	return nu;
}

double MainDistribution::get_mu() const { // получение сдвига
	return mu;
}

double MainDistribution::get_lambda() const { // получение масштаба
	return lambda;
}

void MainDistribution::set_nu(double val) { // установка парметра формы
	nu = val;
}

void MainDistribution::set_mu(double val) { // установка сдвига
	mu = val;
}

void MainDistribution::set_lambda(double val) { // установка масштаба
	lambda = val;
}

double MainDistribution::rand_value() const{
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

double MainDistribution::density(const double& x) const{ // плотность
	if (nu < 0)
	{
		throw 0;
	}
	return 1 / (lambda * pow(2, nu - 1) * (gamma(nu / 2.) * gamma(nu / 2.) / gamma(nu)) * pow(cosh((x - mu) / lambda), nu));
}

double MainDistribution::mathematicalExpectation() const{
   if (nu < 0)
      throw;
   return mu;
}

double MainDistribution::dispersion() const {
   if (nu < 0)
      throw;
   return 0.5 * trigamma(nu/2) * pow(lambda, 2);   
}

double MainDistribution::asymmetry() const {
   if (nu < 0)
      throw;
   return 0;
}

double MainDistribution::kurtosis() const {
   if (nu < 0)
      throw;
   return (0.5 * pentagamma(nu/2))/pow(trigamma(nu/2), 2);
}

void MainDistribution::plot(){
   int N = 1000;
   double x1 = -10;
   double x2 = 10;
   double dx = (x2 - x1) / (double)N;
   string filename = "../../res/distribution_fplot";
   filename.append(to_string(n)).append("_").append(to_string(numerator)).append(".txt");
   ofstream file(filename.c_str());
   for (int i = 0; i < N; i++) {

      file << x1 + i * dx << '\t' << MainDistribution::density(x1 + i * dx) << endl;
   }
   numerator++;
   file.close();
}

string MainDistribution::save_params() const { // сохранение атрибутов в файл
	string filename = "../../distributions/distribution";
	filename.append(to_string(n)).append(".txt");
	ofstream file(filename.c_str());
	file << nu << ' ' << mu << ' ' << lambda << ' ' << n << endl;
	file.close();
	return filename;
}

void MainDistribution::load_params(const std::string& filename) { // �������� ��������� �� �����
	std::ifstream file(filename.c_str());
	file >> nu >> mu >> lambda >> n >> numerator;
	file.close();
}