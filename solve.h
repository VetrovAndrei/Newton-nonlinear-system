#include "func.h"
#pragma once
class solve
{
private: 
	func Fun;
	double B, eps1, eps2, normF0, normFk, normFk1;
	int maxiter;
	std::vector<double> x;
	std::vector<double> F;
	std::vector<double> dx;
	std::vector<double> xk;
	std::vector<std::vector<double>> A;
public:
	void gauss();
	double norm(std::vector<double> X);
	void jacoby(bool dif);
	void sortMatrix(int k);
	void sumElements(int k);
	int iterationB();
	void snu(bool method, bool dif);
	void print(int iter1, std::ofstream &out);
	void print(std::ofstream &out);
	solve(std::ifstream &size, std::ifstream &x0);
	~solve(void);
};

