#include "emp_distribution.hpp"

using namespace std;

// конструкторы
EmpDistribution::EmpDistribution(int N, double* input_samples, int Num) : num(Num), numerator(0) {
	set_samples(N, input_samples);
}

EmpDistribution::EmpDistribution(const std::vector<double>& input_samples, int Num) : num(Num), numerator(0) {
	set_samples(input_samples);
}

EmpDistribution::EmpDistribution(int N, const IDistribution& M, int Num) : num(Num), numerator(0) {
	set_samples(N, M);
}

EmpDistribution::EmpDistribution(FILE* file) { // конструктор из FILE
	if (file != NULL) {
		fscanf_s(file, "%d", &n);
		samples.resize(n);
		for (int i = 0; i < n; i++) {
			fscanf_s(file, "%lf", &samples[i]);
		}
	}
	set_density_segments();
}

EmpDistribution::EmpDistribution(std::string path) { // конструктор из пути std:: string
	FILE* file;
	fopen_s(&file, path.c_str(), "r");
	if (file != NULL) {
		fscanf_s(file, "%d", &n);
		samples.resize(n);
		for (int i = 0; i < n; i++) {
			fscanf_s(file, "%lf", &samples[i]);
		}
	}
	set_density_segments();
}

EmpDistribution::EmpDistribution(const char* path) { // конструктор из пути Си-строки
	FILE* file;
	fopen_s(&file, path, "r");
	if (file != NULL) {
		fscanf_s(file, "%d", &n);
		samples.resize(n);
		for (int i = 0; i < n; i++) {
			fscanf_s(file, "%lf", &samples[i]);
		}
	}
	set_density_segments();
}

EmpDistribution::EmpDistribution(std::ifstream& file) { // конструктор из ifstream
	file >> n;
	samples.resize(n);
	for (int i = 0; i < n; i++) {
		file >> samples[i];
	}
	set_density_segments();
}

EmpDistribution::EmpDistribution(const EmpDistribution& ED) {
	*this = ED;
}

void EmpDistribution::set_samples(const std::vector<double>& input_samples) { // выборка размера n для основного распределения
	n = input_samples.size();
	samples.resize(n);
	for (int i = 0; i < input_samples.size(); i++) {
		samples[i] = input_samples[i];
	}
	set_density_segments();
}

void EmpDistribution::set_samples(int N, double* input_samples) { // выборка размера n для основного распределения
	n = N;
	samples.resize(n);
	for (int i = 0; i < N; i++) {
		samples[i] = input_samples[i];
	}
	set_density_segments();
}

void EmpDistribution::set_samples(int N, const IDistribution& D) { // выборка размера n для 
	if (N > 0) {
		n = N;
		samples.resize(n);
		for (int i = 0; i < N; i++) {
			samples[i] = D.rand_value();
		}
		set_density_segments();
	}
	else throw 0;
}

std::vector<double> EmpDistribution::get_samples() const {
	return samples;
}

std::vector<segment> EmpDistribution::get_density_segments() const {
	return density_segments;
}

double EmpDistribution::rand_value() const {
	if (samples.size() > 0) {
		double r1 = 0;
		double r2 = 0;

		double segment_len = density_segments[0].end - density_segments[0].begin;

		double p = 0;
		double result = 0;

		for (int j = 0; j < n; j++) {
			do r1 = static_cast<double>(rand()) / RAND_MAX; while (r1 == 0 || r1 == 1);
			do r2 = static_cast<double>(rand()) / RAND_MAX; while (r2 == 0);
			p = 0;
			for (int i = 0; i < k; i++) {
				if (p < r1 && r1 < p + density_segments[i].val * segment_len) {
					return density_segments[0].begin + segment_len * (i + r2);
				}
				else p += density_segments[i].val * segment_len;
			}
		}
	}
	else throw;
}

void EmpDistribution::set_density_segments() {
	if (samples.size() == 0) {
		throw 0;
	}

	k = (int)log2(samples.size()) + 1;
	density_segments.resize(k);
	sort(samples.begin(), samples.end());

	double samples_max = samples[samples.size() - 1];
	double samples_min = samples[0];
	double samples_len = samples_max - samples_min;

	for (int i = 0, j = 0, ni = 0; i < k; i++) {
		for (ni = 0; j < samples.size(); j++) { // подсчёт количества элементов в промежутке Xi
			if (samples[j] < samples_min + ((double)(i + 1) / (double)k) * samples_len) { // проверка на принадлежность очередного элемента промежутку Xi
				ni++;
			}
			else {
				break;
			}
		}
		density_segments[i] = { samples_min + ((double)i / (double)k) * samples_len,
			samples_min + ((double)(i + 1) / (double)k) * samples_len,
			(double)ni / (samples.size() * samples_len / (double)k) }; // подсчёт окончен и плотность найдена
	}
}

double EmpDistribution::density(const double& x) const { // samples - выборка из заданного распределения для случайной величины
	if (samples.size() == 0) {
		throw 0;
	}
	
	for (int i = 0; i < k; i++) {
		if (density_segments[i].begin <= x && x < density_segments[i].end) {
			return density_segments[i].val;
		}
	}
	if (abs(x - density_segments[k - 1].end) < 1e-6) return density_segments[k - 1].val;
	else return 0;
}

double EmpDistribution::mathematicalExpectation() const {
	double sum = 0;
	for (int i = 0; i < samples.size(); i++) {
		sum += samples[i];
	}
	return sum / (double)samples.size();
}

double EmpDistribution::dispersion() const {
	double sum = 0;
	double exp_val = mathematicalExpectation();
	for (int i = 0; i < samples.size(); i++) {
		sum += (samples[i] - exp_val) * (samples[i] - exp_val);
	}
	return sum / (double)samples.size();
}

double EmpDistribution::asymmetry() const {
	double sum = 0;
	double exp_val = mathematicalExpectation();
	for (int i = 0; i < samples.size(); i++) {
		sum += (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val);
	}
	return sum / ((double)samples.size() * pow(dispersion(), 1.5));
}

double EmpDistribution::kurtosis() const {
	double sum = 0;
	double exp_val = mathematicalExpectation();
	for (int i = 0; i < samples.size(); i++) {
		sum += (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val);
	}
	return sum / ((double)samples.size() * pow(dispersion(), 2));
}

void EmpDistribution::plot() {
	if (samples.size() > 0) {
		int N = 1000;
		double x1 = -10;
		double x2 = 10;
		double dx = (x2 - x1) / (double)N;
		std::string filename = "../../res/emp_distribution_fplot";
		filename.append(std::to_string(num)).append("_").append(std::to_string(numerator)).append(".txt");
		std::ofstream file(filename.c_str());

		double samples_max = density_segments[k-1].end;
		double samples_min = density_segments[0].begin;

		int ni = 0;
		if (x1 > samples_min) {
			while (ni < k) {
				if (x1 < samples_min + ((double)(ni + 1) / (double)k) * (samples_max - samples_min)) {
					break;
				}
				else ni++;
			}
		}

		for (int i = 0; i < N && ni < k; i++) {
			if (x1 + i * dx < samples_min || x1 + i * dx > samples_max) {
				file << x1 + i * dx << '\t' << 0 << std::endl;
			}
			else if (x1 + i * dx < density_segments[ni].end) {
				file << x1 + i * dx << '\t' << density_segments[ni].val << std::endl;
			}
			else {
				ni++;
				file << x1 + i * dx << '\t' << density_segments[ni].val << std::endl;
			}
		}
		file.close();
		numerator++;
	}
	else throw 0;
}

//сохранение атрибутов в файл
std::string EmpDistribution::save_params() const {
	std::string filename = "../../distributions/emp_distribution";
	filename.append(std::to_string(num)).append(".txt");
	std::ofstream file(filename.c_str());
	file << n << ' ';
	for (int i = 0; i < n; i++)
		file << samples[i] << ' ';
	file.close();
	return filename;
}

void EmpDistribution::load_params(const std::string& filename) {
	std::ifstream file(filename.c_str());
	file >> n;
	samples.resize(n);
	for (int i = 0; i < n; i++)
		file >> samples[i];
	file.close();
	set_density_segments();
}
