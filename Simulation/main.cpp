#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <math.h>
#include <algorithm>
#include <ctime>

std::vector<double> MatrixVecMult(std::vector<std::vector<double> > &, std::vector<double> &);
std::vector<double> RunMC(int);
double L2Norm(std::vector<double> &);



int main()
{
	const int n2steps = 10;
	const int M = 10000;
	

	//Run Calcs and only save end values
	std::vector<double> Xs(n2steps);



	std::cout << "Running in Serial" << std::endl;
	time_t start = time(0);
	for (int m = 0; m < M; m++)
	{
		std::vector<double> X = RunMC(n2steps);
		for (int i = 0; i < n2steps; i++)
			Xs[i] += X[i];
	}
	time_t end = time(0);
	double time = difftime(end, start);
	
	for (int i = 0; i < Xs.size(); i++)
	{
		Xs[i] = std::sqrt(Xs[i] / M);
		std::cout << Xs[i] << ", ";
	}
	std::cout << std::endl;
	std::cout << "Elapsed time in seconds: " << time << "s" << std::endl;
	// std::cin.get();
	return 1;
}


std::vector<double> RunMC(int n2steps_)
{
	// SDE Params
	const double lambda = 2.5;
	const double mu = 1;
	const double T = 1;
	int Nmax;
	std::vector<std::vector<double> > eta = { { 2 / std::sqrt(10), 1 / std::sqrt(10) },{ 1 / std::sqrt(10), 2 / std::sqrt(10) } };
	std::vector<int> Narr(n2steps_);
	for (int i = 0; i < n2steps_; i++) { Narr[i] = (int)std::pow(2, i); }
	Nmax = Narr.back();

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.0);
	std::vector<double> Xs(n2steps_);
	std::vector<std::vector<double> > Xn(n2steps_, std::vector<double>(2,1)); // Initialise matrix of ones
	
	//Generate Random Number Arrays
	std::vector<double> randn1(Nmax);
	std::vector<double> randn2(Nmax);
	for (int i = 0; i < Nmax; i++)
	{
		randn1[i] = distribution(generator);
		randn2[i] = distribution(generator);
	}

	for (int n = 0; n < Nmax; n++)
	{
		for (int i = 0; i < n2steps_; i++)
		{
			if (n % Narr[n2steps_ - 1 - i] == 0)
			{
				std::vector<double> X = Xn.at(i);
				std::vector<double> z = { randn1[n], randn2[n] };
				double h = T / Narr[i];
				double tamedCoeff = 1 / (1 + std::pow(n, -1 / 2) * L2Norm(X));
				std::vector<double> etaZ = MatrixVecMult(eta, z);
				double l2nX = L2Norm(X);

				for (int j = 0; j < X.size(); j++)
				{
					Xn[i][j] += X[j] * (lambda * (mu - l2nX)) * h * tamedCoeff + etaZ[j] * tamedCoeff * std::pow(l2nX, 3 / 2) * std::sqrt(h);

				}
			}
		}
	}

	std::vector<double> XnLastCol = Xn.back();

	for (int i = 0; i < n2steps_; i++)
	{
		std::vector<double> tmp(n2steps_);
		
		std::transform(XnLastCol.begin(), XnLastCol.end(), Xn[i].begin(), tmp.begin(), std::minus<double>());
		Xs[i] = std::pow(L2Norm(tmp), 2);
	}

	return Xs;
}

double L2Norm(std::vector<double> & arr)
{
	double sum = 0;
	for (int i = 0; i < arr.size(); i++)
		sum += arr[i] * arr[i];
	return std::sqrt(sum);
}

std::vector<double> MatrixVecMult(std::vector<std::vector<double> > &Matrix, std::vector<double> &Vec)
{

	int mRows = Matrix.size();
	int mCols = Matrix[0].size();
	int vSize = Vec.size();
	std::vector<double> result(vSize);

	for (int r = 0; r < mRows; r++)
	{
		for (int c = 0; c < mCols; c++)
		{
			double tmp = 0;
			for (int v = 0; v < vSize; v++)
			{
				tmp += Matrix[r][v] * Vec[v];
			}
			result[r] = tmp;
		}
	}

	return result;
}