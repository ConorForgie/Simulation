#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <omp.h>
#include <chrono>


std::vector<double> MatrixVecMult(std::vector<std::vector<double> > &, std::vector<double> &);
std::vector<double> RunMC(int, int);
double L2Norm(std::vector<double> &);

// To use the divide and conquer paradigm - we run many simulations at once, save the resultant vectors and then combine the results into one final vector.
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = omp_orig)


int main()
{
	// Testing
	/*std::vector<std::vector<double> > mat = { {2,0},{0,2} };
	std::vector<double> v = { 3,4 };
	std::vector<double> res = MatrixVecMult(mat, v);
	for (int i = 0; i < res.size(); i++)
		std::cout << res[i] << std::endl;
	double l2n = L2Norm(v);
	std::cout << "L2norm " << l2n << std::endl;*/

	// stepsize = 2^-n2start decreasing to 2^-n2steps
	const int n2start = 6;
	const int n2steps = 20;
	// Run n_trials simulations
	const int n_trials = 1000;
	

	//Run Calcs and only save end values - we don't care about the history of the path
	std::vector<double> Xs(n2steps-n2start);

	// Check if running in parallel or serial
	#ifdef _OPENMP
		std::cout << "Running in parallel" << std::endl;
    #else
		std::cout << "Running in Serial" << std::endl;
	#endif


	std::cout << "Number of simulations = " << n_trials << std::endl;
	std::cout << "Maximum Number of Steps = 2^" << n2steps-1 << std::endl;

	time_t start = time(0);

	int num_threads = 0;

	#ifdef _OPENMP
		double wtime = omp_get_wtime();
	#endif

	// The main body of work - The vector Xs will keep the results and is the target vector for the reduction
	#pragma omp parallel for reduction(vec_double_plus : Xs)
	for (int m = 0; m < n_trials; m++)
	{
		#ifdef _OPENMP
		if(num_threads == 0)
			num_threads = omp_get_num_threads();
		#else
			num_threads = 1;
		#endif
		// Run the Monte Carlo Simulation and collect the results for each stepsize
		std::vector<double> X = RunMC(n2start,n2steps);
		// Collect the results of each simulation
		for (int i = 0; i < n2steps - n2start; i++)
			Xs[i] += X[i];
	}
	
	#ifdef _OPENMP
		// Calculate the wall (run) time
		double fwtime = omp_get_wtime() - wtime;
	#endif

	std::cout << "Number of threads = " << num_threads << std::endl << "---------------------------------" << std::endl << std::endl;

	time_t end = time(0);
	double time = difftime(end, start);
	
	// Build the array of number of steps i.e. 2^n2start to 2^n2steps
	std::vector<double> Narr(n2steps - n2start);
	for (int i = 0; i < n2steps - n2start; i++) { Narr[i] = (int)std::pow(2, i + n2start); }

    // Console Output - the error results
	for (int i = 0; i < Xs.size(); i++)
	{
		// Need to average the results vector Xs
		Xs[i] = std::sqrt(Xs[i] / n_trials);
		std::cout << "2^"<< i+n2start << " = " << Xs[i] << std::endl;
	}
	std::cout << std::endl;
	#ifdef _OPENMP
		std::cout << "Elapsed wtime in seconds: " << fwtime << "s" << std::endl;	
	#endif
	std::cout << "Elapsed time_t in seconds: " << time << "s" << std::endl;
	//std::cin.get(); // For VS
	return 1;
}

// The Monte Carlo Engine
std::vector<double> RunMC(int n2start_, int n2steps_)
{
	// SDE Params
	const double lambda = 2.5;
	const double mu = 1;
	const double T = 1;

	int Nmax;
	int n2length = n2steps_ - n2start_ + 1;
	std::vector<std::vector<double> > eta = { { 2 / std::sqrt(10), 1 / std::sqrt(10) },{ 1 / std::sqrt(10), 2 / std::sqrt(10) } };
	
	std::vector<int> Narr(n2length);
	for (int i = 0; i < n2length; i++) { Narr[i] = (int)std::pow(2, i + n2start_); }
	Nmax = Narr.back();

	std::default_random_engine generator;
	generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::normal_distribution<double> distribution(0.0, 1.0);
	std::vector<std::vector<double> > Xn(n2length, std::vector<double>(2,1)); // Initialise matrix of ones
	
	//Generate Random Number Arrays
	std::vector<double> randn1(Nmax);
	std::vector<double> randn2(Nmax);
	for (int i = 0; i < Nmax; i++)
	{
		randn1[i] = std::sqrt(T / Nmax) * distribution(generator);
		randn2[i] = std::sqrt(T / Nmax) * distribution(generator);
	}

	for (int n = 0; n < Nmax; n++)
	{
		for (int i = 0; i < n2length; i++)
		{
			if (n % Narr[n2length - 1 - i] == 0)
			{
				double r1 = 0; double r2=0;
				for(int j =0; j< Narr[n2length - 1 - i]; j++)
				{
					r1 += randn1[n+j];
					r2 += randn2[n+j];
				}
				std::vector<double> z = { r1, r2 };
				double h = T / Narr[i];
				double l2nX = L2Norm(Xn[i]);
				double tamedCoeff = 1 / (1 + std::pow(n, -1 / 2) * l2nX);
				std::vector<double> etaZ = MatrixVecMult(eta, z);

				for (int j = 0; j < Xn[i].size(); j++)
				{
					Xn[i][j] += Xn[i][j] * (lambda * (mu - l2nX)) * h * tamedCoeff + etaZ[j] * tamedCoeff * std::pow(l2nX, 3 / 2);

				}
			}
		}
	}

	std::vector<double> X_fine = Xn.back();
	std::vector<double> X_soln(n2length-1);

	for (int i = 0; i < n2length-1; i++)
	{
		for (int j = 0; j < Xn[0].size(); j++)
			Xn[i][j] = X_fine[j] - Xn[i][j];

		//std::transform(XnEnd.begin(), XnEnd.end(), Xn[i].begin(), tmp.begin(), std::minus<double>());
		X_soln[i] = std::pow(L2Norm(Xn[i]),2);
	}
	
	// return the square of the L2norm error
	return X_soln;
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