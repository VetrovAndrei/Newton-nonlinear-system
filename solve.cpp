#include "solve.h"


solve::solve(std::ifstream &size, std::ifstream &x0)
{
	Fun.read(size);
	size >> eps1 >> eps2 >> maxiter;
	x.resize(Fun.n);
	F.resize(Fun.m);
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
	Fun.FullF(F,x);
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

void solve::jacoby(bool dif)
{
	for (int i = 0; i < Fun.m; i++)
	{
		for (int j = 0; j < Fun.n; j++)
		{
			if (dif)
				A[i][j] = Fun.chisldF(i, j, x);
			else 
				A[i][j] = Fun.dF(i, j, x);
		}
	}
}

void solve::sortMatrix(int k)
{
	double min;
	int imin;
	for (int i = 1; i <= k; i++)
	{
		imin = 0;
		min = F[0];
		for (int j = 0; j < Fun.m - i; j++)
		{
			if ( F[j] < min)
			{
				min = F[j];
				imin = j;
			}
		}
		A[imin].swap(A[Fun.m - i]);
		std::swap(F[imin],F[Fun.m - i]);
	}
}

void solve::sumElements(int k)
{
	int i = Fun.m - k;
	for (int j = 0; j < Fun.n; j++)
	{
		A[i][j] *= 2.0;
	}
	F[i] = pow(F[i],2.0);
	for (int k = i + 1; k < Fun.m; k++)
	{
		for (int j = 0; j < Fun.n; j++)
		{
			A[i][j] += 2.0 * A[k][j];
		}
		F[i] += pow(F[k],2.0);
	}
}

int solve::iterationB()
{
	B = 1;
	bool exit = 1;
	normFk = normFk1;
	for (int i = 0; i < maxiter && exit != 0; i++)
	{
		if(normFk1 < normFk && B < eps1)
		{
			exit = 0;
			break;
		}
		for (int j = 0; j < Fun.n; j++)
		{
			xk[j] = x[j] + dx[j] * B;
		}
		Fun.FullF(F,xk);
		normFk1 = norm(F);
		B /= 2.0;
	}
	normFk = normFk1;
	x = xk;
}

void solve::snu(bool method, bool dif)
{
	std::ofstream out("output.txt");
	int k = Fun.m - Fun.n;
	normFk = normF0;
	int exit = 1;
	int iterB;
	for (int i = 0; i < maxiter && exit != 0; i++)
	{

		jacoby(dif);
		sortMatrix(k + method);
		if (method)
		{
			sumElements(k + method);
		}
		gauss();
		iterB = iterationB();
	}

}

void solve::print(int iter1, int iter2, std::ofstream &out)
{
	out << iter1 << "\t" << iter2 << "\t" << B << "\t" << normFk / normF0 << std::endl;
	for (int i = 0; i < Fun.m; i++)
	{
		out << x[i] << "\t";
	}
	out << std::endl;
}