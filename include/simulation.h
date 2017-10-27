#include <iostream>
#include <cmath>
#include <random>


#ifndef SIMULATION_H
#define SIMULATION_H

class simulation
{
	private:
        // SDE Params
        const double lambda, mu, T;
        double eta[2][2];
        const int n2steps;
        int Narr[n2steps];
        int Nmax;
        const int M = 1000;
		
	public:
		simulation();
		
        static double [] RunMC();
		static double L2Norm(double [] );
        static double [] MatrixVecMult(double [,] ,double [] );
	
	
}


#endif