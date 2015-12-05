#include "solve.h"


solve::solve(std::ifstream &size, std::ifstream &x0)
{
	Fun.read(size);
	size >> eps1 >> eps2 >> maxiter;
	x.resize(Fun.n);
	F.resize(Fun.n);
	dx.resize(Fun.n);
	xk.resize(Fun.n);
	A.resize(Fun.m);
	for (int i = 0; i < Fun.m; i++)
	{
		A[i].resize(Fun.n);
	}
	for (int i = 0; i < Fun.n; i++)
	{
		x0 >> x[i];
	}
	Fun.F1(F,x);
	normF0 = norm(F);
}


solve::~solve(void)
{
}

double solve::norm(std::vector<double> X)
{
	double sum = 0;
	for (int i = 0; i < X.size(); i++)
	{
		sum += pow(X[i],2.0);
	}
	sum = sqrt(sum);
	return sum;
}

void solve::gauss()
{
	for (int i = 0; i < Fun.n; i++)
	{
		dx[i] = -F[i];
	}
    for (int i = 0; i<Fun.n; i++)
    {
		double max = A[i][i];
        int k = i;
        for (int c = i; c < Fun.n; c++)
        {
			if (A[c][i] > max)
            {
				max = A[c][i];
                k = c;
            }
         }
         A[i].swap(A[k]);
		 std::swap(dx[i],dx[k]);
         for (int j=i+1; j<Fun.n; j++)
         {
			double m = A[j][i]/A[i][i];
            for (int k = 0; k < Fun.n; k++)
			{
				A[j][k]=A[j][k]-m*A[i][k];
			}
            dx[j] = dx[j] - m * dx[i]; 
         }
    }
    for (int k = Fun.n-1; k >= 0; k--)
    {
		double buf = 0;
        for (int j = k+1; j < Fun.n; j++)
        {
			buf += A[k][j]*F[j];
        }
        dx[k] = dx[k] - buf;
        dx[k] = dx[k]/A[k][k];
     }
}

void 