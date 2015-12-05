#include "func.h"
#pragma once
class solve
{
private: 
	func Fun;
	double B, eps1, eps2, normF0, normFk;
	int maxiter;
	std::vector<double> x;
	std::vector<double> F;
	std::vector<double> dx;
	std::vector<double> xk;
	std::vector<std::vector<double>> A;
public:
	void gauss();
	double norm(std::vector<double> X);
	void exclude();
	solve(std::ifstream &size, std::ifstream &x0);
	~solve(void);
};

