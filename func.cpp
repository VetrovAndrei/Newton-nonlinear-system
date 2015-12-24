#include "func.h"


func::func(void)
{
	h = 0.5;
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
	for (int i = 0; i < m; i++)
	{
		F[i] = this->F(i, x);
	}
}

double func::F(int i, std::vector<double> &x)
{
	double res = 0;
	switch(i)
	{
	case 0:
		{
			res = pow(x[0]-2,2) + pow(x[1]-2,2) - 4;
			break;
		}
	case 1:
		{
			res = pow(x[0]-5,2) + pow(x[1]-2,2) - 4;
			break;
		}

	}
	return res;
}

double func::dF(int i, int j, std::vector<double> &x)
{
	double res = 0;
	switch(i)
	{
	case 0:
		{
			switch(j)
			{
				case 0:
				{
					res = 2 * pow(x[0]-2,2);
					break;
				}
				case 1:
				{
					res = 2 * pow(x[1]-2,2);
					break;
				}
			}
			break;
		}
	case 1:
		{
			switch(j)
			{
				case 0:
				{
					res = 2 * pow(x[0]-5,2);
					break;
				}
				case 1:
				{
					res = 2 * pow(x[1]-2,2);
					break;
				}
			}
			break;
		}
	}
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