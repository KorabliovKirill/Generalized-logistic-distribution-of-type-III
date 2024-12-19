#pragma once

#include "distribution.hpp"

template<class Distribution1, class Distribution2>
class MixDistribution : public IDistribution, public IPersistent {
	private:
		Distribution1 Main1; // композиция - строгая зависимость
		Distribution2 Main2;
		int numerator;
		int n;
		double p;

	public:
		MixDistribution() :
			Main1(), Main2(), p(0.5), n(0), numerator(0)
		{}
		MixDistribution(const Distribution1& M1, const Distribution2& M2, double P, int N) :
			Main1(M1), Main2(M2), p(P), n(N), numerator(0)
		{}

		Distribution1& get_M1() { return Main1; }
		Distribution2& get_M2() { return Main2; }

		void set_M1(const Distribution1& M1) { Main1 = M1; }
		void set_M2(const Distribution2& M2) { Main2 = M2; }

		void set_p(const double& val) { p = val; }
		double get_p() const { return p; }
		int get_number() const { return n; } // получение номера распределения

		double rand_value() const override {
			if (p < 0 || p > 1) {
				throw 0;
			}

			double r;
			do r = (double)rand() / RAND_MAX; while (r == 0.0 || r == 1.0);

			if (p < r) return Main1.get_rand_value();
			else return Main2.get_rand_value();
		}

		double density(const double& x) const override {
			if (p < 0 || p > 1) {
				throw 0;
			}

			return (1 - p) * Main1.get_density(x) + p * Main2.get_density(x);
		}

		double mathematicalExpectation() const override {
			if (p <= 0 && p >= 1) {
				throw 0;
			}
			return (Main1.get_expected_value() * (1 - p) + Main2.get_expected_value() * p);
		}

		double dispersion() const override {
			if (p <= 0 && p >= 1) {
				throw 0;
			}

			double D1 = (1 - p) * (pow(Main1.get_expected_value(), 2) + Main1.get_dispersion());
			double D2 = p * (pow(Main2.get_expected_value(), 2) + Main2.get_dispersion());
			return D1 + D2 - pow(get_expected_value(), 2);
		}

		double asymmetry() const override {
			if (p <= 0 && p >= 1) {
				throw 0;
			}

			double M1 = Main1.get_expected_value();
			double M2 = Main2.get_expected_value();
			double M = get_expected_value();
			double D1 = Main1.get_dispersion();
			double D2 = Main2.get_dispersion();
			double el1 = (1 - p) * (pow(M1 - M, 3) + 3 * (M1 - M) * D1 + pow(D1, 3.0 / 2.0) * Main1.get_skewness());
			double el2 = p * (pow(M2 - M, 3) + 3 * (M2 - M) * D2 + pow(D2, 3.0 / 2.0) * Main2.get_skewness());
			return (el1 + el2) / pow(get_dispersion(), 2);
		}

		double kurtosis() const override {
			if (p <= 0 && p >= 1) {
				throw 0;
			}

			double M1 = Main1.get_expected_value();
			double M2 = Main2.get_expected_value();
			double M = get_expected_value();
			double D1 = Main1.get_dispersion();
			double D2 = Main2.get_dispersion();
			double el1 = (1 - p) * (pow(M1 - M, 4) + 6 * pow(M1 - M, 2) * D1 + 4 * (M1 - M) * pow(D1, 3.0 / 2.0) * Main1.get_skewness() + pow(D1, 2) * (Main1.get_kurtosis() + 3));
			double el2 = p * (pow(M2 - M, 4) + 6 * pow(M2 - M, 2) * D2 + 4 * (M2 - M) * pow(D2, 3.0 / 2.0) * Main2.get_skewness() + pow(D2, 2) * (Main2.get_kurtosis() + 3));

			return (el1 + el2) / pow(get_dispersion(), 2) - 3;
		}

		void plot() override {
			int N = 1000;
			double x1 = -10;
			double x2 = 10;
			double dx = (x2 - x1) / (double)N;
			std::string filename = "../../res/mix_distribution_fplot";
			filename.append(std::to_string(n)).append("_").append(std::to_string(Main1.get_number())).append("_");
			filename.append(std::to_string(Main2.get_number())).append("_").append(std::to_string(numerator)).append(".txt");
			std::ofstream file(filename.c_str());
			for (int i = 0; i < N; i++) {

				file << x1 + i * dx << '\t' << get_density(x1 + i * dx) << std::endl;
			}
			file.close();
			numerator++;
		}

		std::string save_params() const override { // сохранение атрибутов в файл
			std::string filename = "../../distributions/mix_distribution";
			filename.append(std::to_string(n)).append(".txt");;
			std::ofstream file(filename.c_str());
			std::string component = "";
			component = Main1.save_params();
			file << component;
			component = Main2.save_params();
			file << component;
			file.close();
			return filename;
		}

		void load_params(const std::string& filename) override {
			std::ifstream file(filename.c_str());
			std::string component = "";
			file >> component;
			Main1.load_params(component);
			file >> component;
			Main2.load_params(component);
			file.close();
		}
};