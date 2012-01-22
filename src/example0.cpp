// mainLBmpcTP.cpp
// example file to test simple examples 
// date: November 08, 2011
// author: Xiaojing ZHANG

// matrices are imported from binary file created by MATLAB


#include <iostream>		// I-O
#include <Eigen/Dense>	// matrix computation
#include "LBmpcTP.h"	// class template
#include <sys/time.h>

using namespace Eigen;
using namespace std;

int main()
{
	// ------ SPECIFY parameters ------- 
	const int _N = 10;		// MPC horizon
	const int _m = 2;		// #input
	const int _n = 5;		// #states
	const int _nSt = 10;	// # state constraints
	const int _nInp = 4;	// # state constraints
	const int _nF_xTheta = 10;	// # Omega constraints
	const int _pos_omega = 1;	// <= _N
	const char fileName[] = "ConstrParam.bin";
	bool verbose = 1;	// '0' if it should shut up
	int status;		// stores status code
	int iterations = 1;
	
	
	// -------- object instantiation -----------
	// fileName contains name of file with parameters
	// bool verbose: 1 or 0
	LBmpcTP<double, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> myObj(	fileName, verbose);


	// ----------- SPECIFY arguments for step() ---------------
	// those matrices are computed externely
	// -----------  sizes of matrices ---------------
		Matrix<double, _n, _n> Lm;		// n x n
		Matrix<double, _n, _m> Mm;		// n x m
		Matrix<double, _n, 1> tm;		// n x 1
		Matrix<double, _n, 1> x_hat;		// n x 1, state estimate
		Matrix<double, _n, 1> x_star[_N];	// n x 1, [_N], tracking
		Matrix<double, _m, 1> u_opt;			// m x 1, optimal input is saved there
	
	// -------- they are updated at each time step ------------
	Lm << 1,  2,  0,  1,  2,
	      -2,  1.3,  1,  2,  2.3,
	       1,  2, -1,  0,  -1,
	       1,  2,  2,  -2.3,  1,
	      0, 0,  2,  1.4, -2;
	
	Mm << 1,  1.4,
	   2,  -1,
	   1,   2,
	   0,   0,
	   2, -1;
	
	tm << -1, 2, 1, 1, 2;
	x_hat << 3, 3, -2, 3, 4;

	for(int i = 0; i <= _N-1; i++)
	{
		x_star[i].setZero();
	}
	
	struct timeval start;
	struct timeval end;
	double elapsedTime = 0;
	double timeTmp;
	int newtonSteps = 0;
	gettimeofday(&start, NULL);	
	for (int i = 0; i < iterations; i++)
	{
		status = myObj.step( Lm, Mm, tm, x_hat, x_star );
		newtonSteps += myObj.n_iter_last;
		if (status)
		{
			cerr << "status error at iteration " << i << ": " << status << endl;
			return 1;
		}
	}
	gettimeofday(&end, NULL);
	timeTmp =   (end.tv_sec*1000.0 + end.tv_usec/1000.0) - (start.tv_sec*1000.0 + start.tv_usec/1000.0);
	elapsedTime += timeTmp;
	
	cout << "elapsedTime with " << iterations << " iterations: " << elapsedTime << " [ms]" << endl;
	cout << "needed " << newtonSteps << " Newton steps" << endl;
	
	return 0;
}



















