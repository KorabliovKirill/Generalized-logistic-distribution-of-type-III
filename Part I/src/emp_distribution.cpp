#include "emp_distribution.hpp"

namespace emp_distribution{
   vector<double> samples(int n, double nu, double mu, double lambda){
      if (n > 0) {
			vector<double> samples;
			for (int i = 0; i < n; i++) {
				samples.emplace_back(distribution::rand_value(nu, mu, lambda));
			}
			return samples;
		}
		else throw;
   }

	vector<double> samples(int n, double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p){
      if (n > 0) {
			vector<double> samples;
			for (int i = 0; i < n; i++) {
				samples.emplace_back(mix_distribution::rand_value(nu1, nu2, mu1, mu2, lambda1, lambda2, p));
			}
			return samples;
		}
		else throw;
   }

	vector<double> samples(int n, vector<double>& samples){
      if (n > 0) {
			return emp_distribution::rand_values(n, samples);
		}
		else throw;
   }

	vector<double> rand_values(int n, vector<double>& samples){
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
		else throw;
   }

	double density(double x, vector<double>& samples){
      if (samples.size() == 0) {
			return 0;
		}

		double samples_max = *max_element(samples.begin(), samples.end());
		double samples_min = *min_element(samples.begin(), samples.end());

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

	double mathematicalExpectation(vector<double>& samples){
      double sum = 0;
		for (int i = 0; i < samples.size(); i++) {
			sum += samples[i];
		}
		return sum / (double)samples.size();
   }

	double dispersion(vector<double>& samples){
      double sum = 0;
		double exp_val = mathematicalExpectation(samples);
		for (int i = 0; i < samples.size(); i++) {
			sum += (samples[i] - exp_val) * (samples[i] - exp_val);
		}
		return sum / (double)samples.size();
   }

	double assymetry(vector<double>& samples){
      double sum = 0;
		double exp_val = mathematicalExpectation(samples);
		for (int i = 0; i < samples.size(); i++) {
			sum += (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val);
		}
		return sum / ((double)samples.size() * pow(dispersion(samples), 1.5));
   }

	double kurtosis(vector<double>& samples){
      double sum = 0;
		double exp_val = mathematicalExpectation(samples);
		for (int i = 0; i < samples.size(); i++) {
			sum += (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val) * (samples[i] - exp_val);
		}
		return sum / ((double)samples.size() * pow(dispersion(samples), 2));
   }

	void plot(int& n, vector<double>& samples){
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
		else throw;
	}
}