#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
#pragma once
class func
{
public:
	int n, m;
	func(void);
	void read(std::ifstream &size);
	~func(void);
	double F1(int i, std::vector<double> &x);
	void F1(std::vector<double> &F, std::vector<double> &x);
	double dF1(int i, int j, std::vector<double> &x);
};

