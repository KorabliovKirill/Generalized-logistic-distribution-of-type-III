#include "mix_distribution.hpp"

MixDistribution::MixDistribution(FILE* file) {
	double nu = 0.0, mu = 0.0, lambda = 0.0;
	if (file != NULL) {
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
		Main1.set_nu(nu);
		Main1.set_mu(mu);
		Main1.set_lambda(lambda);
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
		Main2.set_nu(nu);
		Main2.set_mu(mu);
		Main2.set_lambda(lambda);
		fscanf_s(file, "%lf", &p);
	}
}
MixDistribution::MixDistribution(string path) {
	FILE* file;
	double nu = 0, mu = 0, lambda = 0;
	fopen_s(&file, path.c_str(), "r");
	if (file != NULL) {
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
		Main1.set_nu(nu);
		Main1.set_mu(mu);
		Main1.set_lambda(lambda);
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
		Main2.set_nu(nu);
		Main2.set_mu(mu);
		Main2.set_lambda(lambda);
		fscanf_s(file, "%lf", &p);
	}
}
MixDistribution::MixDistribution(ifstream& file) {
	double nu = 0.0, mu = 0.0, lambda = 0.0;
	file >> nu >> mu >> lambda;
	Main1.set_nu(nu);
	Main1.set_mu(mu);
	Main1.set_lambda(lambda);
	file >> nu >> mu >> lambda >> p;
	Main2.set_nu(nu);
	Main2.set_mu(mu);
	Main2.set_lambda(lambda);
}
MixDistribution::MixDistribution(const char* path) {
	FILE* file;
	double nu = 0.0, mu = 0.0, lambda = 0.0;
	fopen_s(&file, path, "r");
	if (file != NULL) {
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
		Main1.set_nu(nu);
		Main1.set_mu(mu);
		Main1.set_lambda(lambda);
		fscanf_s(file, "%lf", &nu);
		fscanf_s(file, "%lf", &mu);
		fscanf_s(file, "%lf", &lambda);
		Main2.set_nu(nu);
		Main2.set_mu(mu);
		Main2.set_lambda(lambda);
		fscanf_s(file, "%lf", &p);
	}
}
void MixDistribution::set_p(const double& val) {
	p = val;
}

double MixDistribution::get_p() const {
	return p;
}

double MixDistribution::rand_value(MainDistribution M1, MainDistribution M2, double p) {

	if (M1.get_nu() < 0 || M2.get_nu() < 0 || p < 0 || p > 1) {
		throw 0;
	}

	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0.0 || r == 1.0);

	if (p < r) return M1.rand_value();
	else return M2.rand_value();
}

double MixDistribution::density(double x, MainDistribution M1, MainDistribution M2, double p) {
	if (M1.get_nu() < 0 || M2.get_nu() < 0 || p < 0 || p > 1) {
		throw 0;
	}

	return (1 - p) * M1.density(x) + p * M2.density(x);
}

double MixDistribution::mathematicalExpectation(MainDistribution M1, MainDistribution M2, double p) {
	if (M1.get_nu() < 0 || M2.get_nu() < 0 || p <= 0 && p >= 1) {
		throw 0;
	}
	return (M1.mathematicalExpectation() * (1 - p) + M2.mathematicalExpectation() * p);
}

double MixDistribution::dispersion(MainDistribution M1, MainDistribution M2, double p){
	if (M1.get_nu() < 0 || M2.get_nu() < 0 || p <= 0 && p >= 1) {
		throw 0;
	}

	double D1 = (1 - p) * (pow(M1.mathematicalExpectation(), 2) + M1.dispersion());
	double D2 = p * (pow(M2.mathematicalExpectation(), 2) + M2.dispersion());
	return D1 + D2 - pow(mathematicalExpectation(M1, M2, p), 2);
}

double MixDistribution::assymetry(MainDistribution M1, MainDistribution M2, double p) {
	if (M1.get_nu() < 0 || M2.get_nu() < 0 || p <= 0 && p >= 1) {
		throw 0;
	}

	double Main1 = M1.mathematicalExpectation();
	double Main2 = M2.mathematicalExpectation();
	double Main = mathematicalExpectation(M1, M2, p);
	double D1 = M1.dispersion();
	double D2 = M2.dispersion();
	double el1 = (1 - p) * (pow(Main1 - Main, 3) + 3 * (Main1 - Main) * D1 + pow(D1, 3.0 / 2.0) * M1.asymmetry());
	double el2 = p * (pow(Main2 - Main, 3) + 3 * (Main2 - Main) * D2 + pow(D2, 3.0 / 2.0) * M2.asymmetry());
	return (el1 + el2) / pow(dispersion(M1, M2, p), 2);
}

double MixDistribution::kurtosis(MainDistribution M1, MainDistribution M2, double p) {
	if (M1.get_nu() < 0 || M2.get_nu() < 0 || p <= 0 && p >= 1) {
		throw 0;
	}

	double Main1 = M1.mathematicalExpectation();
	double Main2 = M2.mathematicalExpectation();
	double Main = mathematicalExpectation(M1, M2, p);
	double D1 = M1.dispersion();
	double D2 = M2.dispersion();
	double el1 = (1 - p) * (pow(Main1 - Main, 4) + 6 * pow(Main1 - Main, 2) * D1 + 4 * (Main1 - Main) * pow(D1, 3.0 / 2.0) * M1.asymmetry() + pow(D1, 2) * (M1.kurtosis() + 3));
	double el2 = p * (pow(Main2 - Main, 4) + 6 * pow(Main2 - Main, 2) * D2 + 4 * (Main2 - Main) * pow(D2, 3.0 / 2.0) * M2.asymmetry() + pow(D2, 2) * (M2.kurtosis() + 3));

	return (el1 + el2) / pow(dispersion(M1, M2, p), 2) - 3;
}

void MixDistribution::plot(int& n, MainDistribution M1, MainDistribution M2, double p) {
	int N = 1000;
	double x1 = -10;
	double x2 = 10;
	double dx = (x2 - x1) / (double)N;
	std::string filename = "../../res/mix_distribution_fplot";
	filename.append(std::to_string(n)).append(".txt");
	std::ofstream file(filename.c_str());
	for (int i = 0; i < N; i++) {

		file << x1 + i * dx << '\t' << density(x1 + i * dx, M1, M2, p) << std::endl;
	}
	file.close();
}

string MixDistribution::save_params() { // сохранение атрибутов в файл
	std::string filename = "distributions/distribution";
	filename.append(std::to_string(1)).append(".txt");
	std::ofstream file(filename.c_str());
	file << get_M1().get_nu() << ' ' << get_M1().get_mu() << ' ' << get_M1().get_lambda() << ' ' << get_M2().get_nu() << ' ' << get_M2().get_mu() << ' ' << get_M2().get_lambda() << ' ' << get_p() << std::endl;
	file.close();
	return filename;
}