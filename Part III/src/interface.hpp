#pragma once
#include <string>
#define interface struct

interface IDistribution {// интерфейс распределения
	virtual int get_number() const = 0;
	virtual double rand_value() const = 0;
	virtual double density(const double& x) const = 0;
	virtual double mathematicalExpectation() const = 0;
	virtual double dispersion() const = 0;
	virtual double asymmetry() const = 0;
	virtual double kurtosis() const = 0;
	virtual void plot() = 0;
};

interface IPersistent {// интерфейс персистентного объекта
	virtual std::string save_params() const = 0;
	virtual void load_params(const std::string& filename) = 0;
};
