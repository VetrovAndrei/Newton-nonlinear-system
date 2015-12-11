#include "func.h"


func::func(void)
{
}


func::~func(void)
{
}

void func::read(std::ifstream &size)
{
	size >> n >> m;
}

void func::FullF(std::vector<double> &F, std::vector<double> &x)
{
	for (int i = 0; i < n; i++)
	{
		F[i] = this->F(i, x);
	}
}

double func::F(int i, std::vector<double> &x)
{
	double res = 0;
	return res;
}

double func::dF(int i, int j, std::vector<double> &x)
{
	double res = 0;
	return res;
}

double func::chisldF(int i, int j, std::vector<double> &x)
{
	double res = 0;
	double buf = x[j];
	x[j] += h;
	res += F(i, x);
	x[j] = buf;
	x[j] -= h;
	res -= F(i, x);
	res /= 2 * h;
	x[j] = buf;
	return res;
}