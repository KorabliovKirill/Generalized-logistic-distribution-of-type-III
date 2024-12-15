#pragma once

#include "distribution.hpp"

namespace mix_distribution{

	class MixDistribution {
		private:
			MainDistribution Main1; // РєРѕРјРїРѕР·РёС†РёСЏ - СЃС‚СЂРѕРіР°СЏ Р·Р°РІРёСЃРёРјРѕСЃС‚СЊ
			MainDistribution Main2;
			double p;

		public:
			MixDistribution() : Main1(), Main2() { p = 0.5; }
			MixDistribution(double nu1, double nu2, double mu1, double mu2, double lambda1, double lambda2, double p) :
				Main1(nu1, mu1, lambda1, p),
				Main2(nu2, mu2, lambda2, p) {}
			MixDistribution(FILE* file);
			MixDistribution(string path);
			MixDistribution(ifstream& file);
			MixDistribution(const char* path);

			MainDistribution& get_M1() { return Main1; }
			MainDistribution& get_M2() { return Main2; }

			void set_M1(const MainDistribution& M1);
			void set_M2(const MainDistribution& M2);

			void set_p(const double& val);
			double get_p() const;

			double rand_value(MainDistribution M1, MainDistribution M2, double p);
			double density(double x, MainDistribution M1, MainDistribution M2, double p);
			double mathematicalExpectation(MainDistribution M1, MainDistribution M2, double p);
			double dispersion(MainDistribution M1, MainDistribution M2, double p);
			double assymetry(MainDistribution M1, MainDistribution M2, double p);
			double kurtosis(MainDistribution M1, MainDistribution M2, double p);
			void plot(int& n, MainDistribution M1, MainDistribution M2, double p);

			string save_params();
	};
}