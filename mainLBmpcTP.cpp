// mainLBmpcTP.cpp
// example file to test simple examples 
// date: November 08, 2011
// author: Xiaojing ZHANG

// matrices are imported from binary file created by MATLAB


#include <iostream>		// I-O
#include <Eigen/Dense>	// matrix computation
#include "LBmpcTP.h"	// class template


using namespace Eigen;
using namespace std;

int main()
{
	// ------ SPECIFY parameters ------- 
	const int _N = 4;		// MPC horizon
	const int _m = 2;		// #input
	const int _n = 5;		// #states
	const int _nSt = 10;	// # state constraints
	const int _nInp = 4;	// # state constraints
	const int _nF_xTheta = 10;	// # Omega constraints
	const int _pos_omega = 4;	// <= _N
	const char fileName[] = "ConstrParam.bin";
	bool verbose = 0;	// '0' if it should shut up
	int status;		// stores status code
	
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
	// x_star_arg[0] << 29.1867 ,  20.3181,   36.6838,   10.5584,  -24.6923;
	// x_star_arg[1] << 401.2845,  262.2318,  -73.8211, -285.3312,  391.4877;
	// x_star_arg[2] << -1.4669*1000  ,  0.0849*1000 ,   0.1898*1000 ,   2.8506*1000 ,  -1.8977*1000;
	// x_star_arg[3] << -0.8542*1000  ,  9.2976*1000   , 7.0920*1000 ,  -1.5346*1000 ,   4.6926*1000;
	
	status = myObj.step( Lm, Mm, tm, x_hat, x_star );
	if (!status)
	{
		u_opt = myObj.u_opt;
		cout << "optimal input:" << endl << u_opt << endl << endl;
	}
	
	return 0;
}



















