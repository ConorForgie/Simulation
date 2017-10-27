#include "../include/simulation.h"

simulation::simulation()
{
	lambda = 2.5;
    mu = 1;
    T = 1;
    eta = { { 2 / Math.Sqrt(10), 1 / Math.Sqrt(10) }, { 1 / Math.Sqrt(10), 2 / Math.Sqrt(10) } };
    n2steps = 10;
    M = 1000;
	for (int i = 0; i < n2steps; i++) { Narr[i] = (int)std::pow(2, i); }
	Nmax = Narr.back();
}

int main()
{
	//Run Calcs and only save end values
	double Xs[n2steps];


	Console.WriteLine("Running in Serial");
	for (int m = 0; m < M; m++)
	{
		double X[] = RunMC();
		for(int i=0; i<sizeof(Xs); i++)
			Xs[i] += X[i];
	}

	for(int i=0; i<sizeof(Xs); i++)
	{
		Xs[i] = std::sqrt(Xs[i]/M);
		std::cout << Xs[i] << ", ";
	}
	cout << std::endl
	return 1;
}

double[] simulation::RunMC()
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0,1.0);	
	double Xs[n2steps];
	double Xn[2,n2steps];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < n2steps; j++)
			Xn[i, j] = 1;

	//Generate Random Number Arrays
	double randn1[Nmax];
	double randn2[Nmax];
	for(int i=0; i<Nmax; i++)
	{
		randn1[i] = distribution(generator);
		randn2[i] = distribution(generator);
	}
	
	for (int n = 0; n < Nmax; n++)
	{
		for (int i = 0; i < n2steps; i++)
		{
			if (n % Narr[sizeof(Narr) - 1 - i] == 0)
			{
				double X[] = { Xn[0, i], Xn[1, i] };
				double z[] = { randn1[n], randn2[n] };
				double h = T / Narr[i];
				double tamedCoeff = 1 / (1 + std::pow(n, -1 / 2) * L2Norm(X));
				double etaZ[] = MatrixVecMult(eta,z);
				double l2nX = L2Norm(X);
				
				for(int j=0; j<sizeof(X); j++)
				{
					Xn[j,i] = X[j]*(lambda * (mu - l2nX)) * h * tamedCoeff + etaZ[j] * tamedCoeff * std::pow(l2nX, 3 / 2) * std::sqrt(h);
					
				}							
			}
		}
	}

	double XnLastCol[] = { Xn[0, n2steps - 1], Xn[1, n2steps - 1] };

	for (int i = 0; i < n2steps; i++)
	{
		double tmp[n2steps];
		for(int j=0; j< n2steps; j++)
			tmp[j] = XnLastCol[j] - Xn[j,i];
		Xs[i] = std::pow(L2Norm(tmp), 2);
	}

	return Xs;
}

double simulation::L2Norm(double arr[][)
{
	double sum = 0;
	for (int i = 0; i < sizeof(arr); i++)
		sum += arr[i] * arr[i];
	return std::sqrt(sum);
}
	
static double[] MatrixVecMult(double Matrix[][] ,double Vec[])
{
	
	int mRows = sizeof(Matrix);
	int mCols = sizeof(Matrix[0]);
	int vSize = sizeof(Vec);
	double result[vsize];
	
	for(int r=0; r<mRows; r++)
	{
		for(int c=0; c<mCols;c++)
		{
			double tmp = 0;
			for(int v=0; v<vSize; v++)
			{
				tmp += Matrix[r][v] * Vec[v]; 
			}
			result[r] = tmp;
		}
	}

	return result;
}

	