/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2013, UC Regents
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the University of California nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

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
	const char fileName[] = "matlab/example0/ConstrParam.bin";
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



















