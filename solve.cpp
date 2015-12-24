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
			 if (abs(A[i][i]) < pow(10.0,-15.0))
				 throw 1;
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
			buf += A[k][j]*dx[j];
        }
        dx[k] = dx[k] - buf;
		if (abs(A[k][k]) < pow(10.0,-15.0))
				 throw 1;
        dx[k] = dx[k]/A[k][k];
     }
}

void solve::jacoby(bool dif)
{
	if (dif)
	{
		for (int i = 0; i < Fun.m; i++)
		{
			for (int j = 0; j < Fun.n; j++)
			{
				A[i][j] = Fun.chisldF(i, j, x);
			}
		}
	}
	else
	{
		for (int i = 0; i < Fun.m; i++)
		{
			for (int j = 0; j < Fun.n; j++)
			{
				A[i][j] = Fun.dF(i, j, x);
			}
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
		min = abs(F[0]);
		for (int j = 0; j < Fun.m - i + 1; j++)
		{
			if (abs(F[j]) < min)
			{
				min = abs(F[j]);
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
		A[i][j] *= 2.0 * F[i];
	}
	F[i] = pow(F[i],2.0);
	for (int k = i + 1; k < Fun.m; k++)
	{
		for (int j = 0; j < Fun.n; j++)
		{
			A[i][j] += 2.0 * A[k][j] * F[k];
		}
		F[i] += pow(F[k],2.0);
	}
}

int solve::iterationB()
{
	B = 1;
	bool exit = 1;
	normFk1 = normFk;
	int i;
	for (i = 0; i < maxiter && exit != 0; i++)
	{
		if(normFk1 < normFk)
		{
			exit = 0;
			break;
		}
		if(B < eps1)
		{
			i = maxiter;
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
	return i;
}

void solve::snu(bool method, bool dif)
{
	double nev = 1;
	std::ofstream out("output.txt");
	out.precision(17);
	std::ofstream o("out.txt");
	o.precision(17);
	print(o);
	int k = Fun.m - Fun.n;
	normFk = normF0;
	int exit = 1;
	int iterB = 0;
	if(abs(normF0) < eps1)
		throw 2;
	for (int i = 0; i < maxiter && exit != 0; i++)
	{
		if(nev < eps2 || iterB == maxiter)
		{
			exit = 0;
			break;
		}
		jacoby(dif);
		sortMatrix(k + method);
		if (method)
		{
			sumElements(k + method);
		}
		gauss();
		iterB = iterationB();
		nev = normFk / normF0;
		print(i, out);
		print(o);
	}

}

void solve::print(std::ofstream &out)
{
	for (int i = 0; i < Fun.n; i++)
	{
		out << x[i] << "\t";
	}
	out << std::endl;
}

void solve::print(int iter1, std::ofstream &out)
{
	out << "iter: " << iter1 + 1;
	out << "; B: "  << 2 * B;
	out << "; nev: " << normFk / normF0 << std::endl;
	for (int i = 0; i < Fun.n; i++)
	{
		out << x[i] << "\t";
	}
	out << std::endl;
	out << std::endl;
}