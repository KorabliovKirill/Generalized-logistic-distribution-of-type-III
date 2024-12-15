#include "emp_distribution.hpp"

using namespace std;

namespace emp_distribution {
	// конструкторы
	EmpDistribution::EmpDistribution() {
		n = 1;
		samples.resize(n);
		samples = {1};
	}

	EmpDistribution::EmpDistribution(int n, MainDistribution M1) {
		this->n = n;
		samples.resize(n);
		std::vector<double> set_samples(int n, MainDistribution M1);
	}

	EmpDistribution::EmpDistribution(int n, MainDistribution M1, MainDistribution M2, double p) {
		this->n = n;
		samples.resize(n);
		std::vector<double> set_samples(int n, MainDistribution M1, MainDistribution M2, double p);
	}

	EmpDistribution::EmpDistribution(int n, vector<double>& val) {
		this->n = n;
		samples.resize(n);
		std::vector<double> set_samples(int n, vector<double>&val);
	}
	EmpDistribution::EmpDistribution(FILE* file) { // конструктор из FILE

		if (file != NULL) {
			fscanf_s(file, "%d", &n);
			samples.resize(n);
			for (int i = 0; i < n; i++) {
				fscanf_s(file, "%lf", &samples[i]);
			}
		}
	}

	EmpDistribution::EmpDistribution(string path) { // конструктор из пути std:: string
		FILE* file;
		fopen_s(&file, path.c_str(), "r");
		if (file != NULL) {
			fscanf_s(file, "%d", &n);
			samples.resize(n);
			for (int i = 0; i < n; i++) {
				fscanf_s(file, "%lf", &samples[i]);
			}
		}
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
	}

	EmpDistribution::EmpDistribution(ifstream& file) { // конструктор из ifstream
		file >> n;
		samples.resize(n);
		for (int i = 0; i < n; i++) {
			file >> samples[i];
		}
	}
	//------------------------------------------------------------------------------------------

	vector<double> EmpDistribution::rand_values(int n, vector<double> samples) {
		if (samples.size() > 0) {
			double r1 = 0;
			double r2 = 0;
			int k = (int)log2(samples.size()) + 1;
			std::vector<double> pi(k);

			sort(samples.begin(), samples.end());
			double samples_max = samples[samples.size() - 1];
			double samples_min = samples[0];
			double samples_len = samples_max - samples_min;

			int ni = 0;
			for (int i = 0, j = 0; i < k; i++) {
				for (ni = 0; j < samples.size(); j++) { // подсчёт количества элементов в промежутке Xi
					if (samples[j] < samples_min + ((double)(i + 1) / (double)k) * samples_len) { // проверка на принадлежность очередного элемента промежутку Xi
						ni++;
					}
					else break;
				}
				pi[i] = (double)ni / samples.size(); // подсчёт окончен и относительная частота найдена
			}

			std::vector<double> result;
			int i_cur = 0;
			double p = 0;
			for (int j = 0; j < n; j++) {
				do r1 = static_cast<double>(rand()) / RAND_MAX; while (r1 == 0 || r1 == 1);
				do r2 = static_cast<double>(rand()) / RAND_MAX; while (r2 == 0);
				i_cur = 0;
				p = 0;
				for (int i = 0; i < k; i++) {
					if (p < r1 && r1 < p + pi[i]) {
						i_cur = i;
						break;
					}
					else p += pi[i];
				}
				result.emplace_back((samples_len / (double)k) * (i_cur + r2) + samples_min);
			}
			return result;
		}
		else throw 0;
	}

	vector<double> EmpDistribution::set_samples(int n, MainDistribution M) { // выборка размера n для основного распределения
		if (n > 0) {
			std::vector<double> val;
			for (int i = 0; i < n; i++) {
				val.emplace_back(M.rand_value());
			}
			return samples = val;
		}
		else throw 0;
	}

	vector<double> EmpDistribution::set_samples(int n, MainDistribution M1, MainDistribution M2, double p) { // выборка размера n для смеси распределений
		if (n > 0) {
			mix_distribution::MixDistribution MD;
			std::vector<double> val;
			for (int i = 0; i < n; i++) {
				val.emplace_back(MD.rand_value(M1, M2, p));
			}
			return samples = val;
		}
		else throw 0;
	}

	vector<double> EmpDistribution::set_samples(int n, vector<double> val) { // выборка размера n для эмпирического распределения
		if (n > 0) {
			return samples = rand_values(n, val);
		}
		else throw 0;
	}
	//--------------------------------------------------------------------
	
	std::vector<double> EmpDistribution::get_samples()
	{
		return samples;
	}

	int EmpDistribution::set_n(int val)
	{
		return n = val;
	}
	int EmpDistribution::get_n(){
		return n;
	}
	//---------------------------------------------------------------------

	double EmpDistribution::density(double x, vector<double> samples) { // samples - выборка из заданного распределения для случайной величины

		if (samples.size() == 0) {
			return 0;
		}

		double samples_max = *std::max_element(samples.begin(), samples.end());
		double samples_min = *std::min_element(samples.begin(), samples.end());

		if (x < samples_min || x > samples_max) return 0;

		double samples_len = samples_max - samples_min;
		int k = (int)log2(samples.size()) + 1;

		sort(samples.begin(), samples.end());

		int ni = 0;
		for (int i = 0; i < k; i++) {
			if (samples_min + ((double)i / (double)k) * samples_len <= x && samples_min + ((double)(i + 1) / (double)k) * samples_len > x) { // проверка на попадание в промежуток Xi
				int j = 0;
				while (samples[j] < samples_min + ((double)i / (double)k) * samples_len) { // поиск номера первого элемента выборки, входящего в Xi
					j++;
				}
				for (ni = 0; ni < samples.size(); j++) { // подсчёт количества элементов в промежутке Xi
					if (samples[j] < samples_min + ((double)(i + 1) / (double)k) * samples_len) { // проверка на принадлежность очередного элемента промежутку Xi
						ni++;
					}
					else break; // если очередной элемент не принадлежит Xi, то подсчёт окончен и можно вернуть результат
				}
				return (double)ni / (samples.size() * samples_len / (double)k);
			}
		}

		return 0;
	}
	double EmpDistribution::mathematicalExpectation(vector<double> samples) {
		double sum = 0;
		for (int i = 0; i < samples.size(); i++) {
			sum += samples[i];
		}
		return sum / (double)samples.size();

	}

	double EmpDistribution::dispersion(vector<double> samples) {
		double sum = 0;
		double exp_val = mathematicalExpectation(samples);
		for (int i = 0; i < samples.size(); i++) {
			sum += (samples[i] - exp_val) * (samples[i] - exp_val);
		}
		return sum / (double)samples.size();
	}

	double EmpDistribution::assymetry(vector<double> samples) {
		double sum = 0;
		double exp_val = mathematicalExpectation(samples);
		for (int i = 0; i < samples.size(); i++) {
			sum += (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val);
		}
		return sum / ((double)samples.size() * pow(dispersion(samples), 1.5));
	}

	double EmpDistribution::kurtosis(vector<double> samples) {
		double sum = 0;
		double exp_val = mathematicalExpectation(samples);
		for (int i = 0; i < samples.size(); i++) {
			sum += (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val);
		}
		return sum / ((double)samples.size() * pow(dispersion(samples), 2));
	}

	void EmpDistribution::plot(int n, std::vector<double> samples) {
		if (samples.size() > 0) {
			int N = 1000;
			double x1 = -10;
			double x2 = 10;
			double dx = (x2 - x1) / (double)N;
			std::string filename = "../../res/emp_distribution_fplot";
			filename.append(std::to_string(n)).append(".txt");
			std::ofstream file(filename.c_str());
			double r1 = 0;
			double r2 = 0;
			int k = (int)log2(samples.size()) + 1;
			std::vector<double> fi(k);

			sort(samples.begin(), samples.end());
			double samples_max = samples[samples.size() - 1];
			double samples_min = samples[0];
			double samples_len = samples_max - samples_min;

			int ni = 0;
			for (int i = 0, j = 0; i < k; i++) {
				for (ni = 0; j < samples.size(); j++) { // подсчёт количества элементов в промежутке Xi
					if (samples[j] < samples_min + ((double)(i + 1) / (double)k) * samples_len) { // проверка на принадлежность очередного элемента промежутку Xi
						ni++;
					}
					else break;
				}
				fi[i] = (double)ni / (samples.size() * samples_len / (double)k); // подсчёт окончен и плотность найдена
			}

			ni = 0;
			if (x1 > samples_min) {
				while (ni < k) {
					if (x1 < samples_min + ((double)(ni + 1) / (double)k) * samples_len) {
						break;
					}
					else ni++;
				}
			}

			for (int i = 0; i < N; i++) {
				if (x1 + i * dx < samples_min || x1 + i * dx > samples_max) {
					file << x1 + i * dx << '\t' << 0 << std::endl;
				}
				else if (x1 + i * dx < samples_min + ((double)(ni + 1) / (double)k) * samples_len) {
					file << x1 + i * dx << '\t' << fi[ni] << std::endl;
				}
				else {
					ni++;
					file << x1 + i * dx << '\t' << fi[ni] << std::endl;
				}
			}
			file.close();
		}
		else throw 0;
	}

	//сохранение атрибутов в файл
	std::string EmpDistribution::save_params()
	{
		std::string filename = "distributions/distribution";
		filename.append(std::to_string(2)).append(".txt");
		std::ofstream file(filename.c_str());
		file << n << std::endl;
		for (int i = 0; i < n; i++)
			file << samples[i]<<' ';
		file.close();
		return filename;
	}
}
