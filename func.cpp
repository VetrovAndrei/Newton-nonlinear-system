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

void func::F1(std::vector<double> &F, std::vector<double> &x)
{
	for (int i = 0; i < n; i++)
	{
		F[i] = F1(i, x);
	}
}

double func::F1(int i, std::vector<double> &x)
{

}

double func::dF1(int i, int j, std::vector<double> &x)
{

}
