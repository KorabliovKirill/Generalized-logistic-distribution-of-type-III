#pragma once

#include "distribution.hpp"
#include "mix_distribution.hpp"
#include <algorithm>

namespace emp_distribution{
	class EmpDistribution{
		private:
			int n; // число элементов выборки
			vector<double> samples; //массив с выборкой

		public:
			EmpDistribution(); // конструктор по умолчанию
			EmpDistribution(int n, MainDistribution M1); // конструктор для стандартного распределения
			EmpDistribution(int n, MainDistribution M1, MainDistribution M2, double p);// конструктор для смеси
			EmpDistribution(int n, vector<double>& val);
			EmpDistribution(FILE* file); // конструктор из FILE
			EmpDistribution(ifstream& file); // конструктор из ifstream
			EmpDistribution(string path); // конструктор из пути std::string
			EmpDistribution(const char* path); // конструктор из пути Си-строки
			
			int set_n(int val);
			int get_n();

			vector<double> set_samples(int n, MainDistribution M1);// создание выборки для стандартного распределения
			vector<double> set_samples(int n, MainDistribution M1, MainDistribution M2, double p);// создание выборки для смеси
			vector<double> set_samples(int n, vector<double> val);// создание выборки для эмпирического распределения
			vector<double> get_samples();// получение выборки

			vector<double> rand_values(int n, vector<double> samples);
			double density(double x, vector<double> samples);
			double mathematicalExpectation(vector<double> samples);
			double dispersion(vector<double> samples);
			double assymetry(vector<double> samples);
			double kurtosis(vector<double> samples);
			void plot(int n, vector<double> samples);

			string save_params();//сохранение атрибутов в файл

			// деструктор, очищающий массив
			~EmpDistribution() {
				vector<double>().swap(samples);
			}
	};
}