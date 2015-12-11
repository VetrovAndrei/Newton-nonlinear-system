#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
#pragma once
class func
{
public:
	int n, m;
	double h;
	func(void);
	void read(std::ifstream &size);
	~func(void);
	double F(int i, std::vector<double> &x);
	void FullF(std::vector<double> &F, std::vector<double> &x);
	double dF(int i, int j, std::vector<double> &x);
	double chisldF(int i, int j, std::vector<double> &x);
};

