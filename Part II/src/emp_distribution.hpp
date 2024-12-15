#pragma once

#include "distribution.hpp"
#include "mix_distribution.hpp"
#include <algorithm>

namespace emp_distribution{
	class EmpDistribution{
		private:
			int n; // ����� ��������� �������
			vector<double> samples; //������ � ��������

		public:
			EmpDistribution(); // ����������� �� ���������
			EmpDistribution(int n, MainDistribution M1); // ����������� ��� ������������ �������������
			EmpDistribution(int n, MainDistribution M1, MainDistribution M2, double p);// ����������� ��� �����
			EmpDistribution(int n, vector<double>& val);
			EmpDistribution(FILE* file); // ����������� �� FILE
			EmpDistribution(ifstream& file); // ����������� �� ifstream
			EmpDistribution(string path); // ����������� �� ���� std::string
			EmpDistribution(const char* path); // ����������� �� ���� ��-������
			
			int set_n(int val);
			int get_n();

			vector<double> set_samples(int n, MainDistribution M1);// �������� ������� ��� ������������ �������������
			vector<double> set_samples(int n, MainDistribution M1, MainDistribution M2, double p);// �������� ������� ��� �����
			vector<double> set_samples(int n, vector<double> val);// �������� ������� ��� ������������� �������������
			vector<double> get_samples();// ��������� �������

			vector<double> rand_values(int n, vector<double> samples);
			double density(double x, vector<double> samples);
			double mathematicalExpectation(vector<double> samples);
			double dispersion(vector<double> samples);
			double assymetry(vector<double> samples);
			double kurtosis(vector<double> samples);
			void plot(int n, vector<double> samples);

			string save_params();//���������� ��������� � ����

			// ����������, ��������� ������
			~EmpDistribution() {
				vector<double>().swap(samples);
			}
	};
}