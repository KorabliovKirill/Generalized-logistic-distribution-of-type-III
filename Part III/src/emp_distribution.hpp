#pragma once

#include "distribution.hpp"
#include "mix_distribution.hpp"
#include <algorithm>

struct segment {
	double begin;
	double end;
	double val;
};

class EmpDistribution : public IDistribution, public IPersistent {
private:
	int n; //����� ��������� �������
	int k; //����� ���������� �������
	std::vector<double> samples; //�������
	std::vector<segment> density_segments;
	int num;
	int numerator;
	void set_density_segments();

public:
	EmpDistribution(): // ����������� �� ���������
		n(1), k(1), samples({ 1 }), density_segments({ {1, 1, 1} }), num(1), numerator(0)
	{}; 
	EmpDistribution(int N, double* input_samples, int Num); // ����������� �� �������
	EmpDistribution(const std::vector<double>& input_samples, int Num);
	EmpDistribution(int N, const IDistribution& D, int Num); // ����������� ��� ������������ �������������
	EmpDistribution(FILE* file); // ����������� �� FILE
	EmpDistribution(std::ifstream& file); // ����������� �� ifstream
	EmpDistribution(std::string path); // ����������� �� ���� std::string
	EmpDistribution(const char* path); // ����������� �� ���� ��-������
	EmpDistribution(const EmpDistribution& emp);

	EmpDistribution& operator=(const EmpDistribution& ED) {
		if (this == &ED) return *this;
		if (n != ED.n) {
			std::vector<double>().swap(samples);
			set_samples(ED.samples);
		}
		if (k != ED.k) {
			std::vector<segment>().swap(density_segments);
			k = ED.k;
			density_segments.resize(k);
			for (int i = 0; i < k; i++) {
				density_segments[i] = ED.density_segments[i];
			}
		}
		num = ED.num;
		numerator = ED.numerator;
		return *this;
	}

	void set_samples(const std::vector<double>& input_samples);// �������� ������� ��� ������������ �������������
	void set_samples(int N, double* input_samples);// �������� ������� ��� ������������ �������������
	void set_samples(int N, const IDistribution& D);// �������� ������� ��� ������������ �������������

	int get_n() const { return n; }// ��������� �������
	std::vector<double> get_samples() const;// ��������� �������
	std::vector<segment> get_density_segments() const;
	int get_number() const { return num; } // ��������� ������ �������������

	double rand_value() const override;
	double density(const double& x) const override;
	double mathematicalExpectation() const override;
	double dispersion() const override;
	double asymmetry() const override;
	double kurtosis() const override;
	void plot() override;

	std::string save_params() const override; //���������� ��������� � ����
	void load_params(const std::string& filename) override; //�������� ��������� � ����

	// ����������, ��������� ������
	~EmpDistribution() {
		std::vector<double>().swap(samples);
		std::vector<segment>().swap(density_segments);
	}
};