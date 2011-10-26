// mainLBmpcTP.cpp
// example file to test simple examples 
// date: October 17, 2011
// author: Xiaojing ZHANG
// 
// horizon: N = 4
// states: n = 5
// input: m = 2;

// matrices are imported from binary file created by MATLAB


#include <iostream>		// I-O
#include <fstream>		// read binary data
#include <Eigen/Dense>	// matrix computation
#include "LBmpcTP.h"	// class template

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
	const int _pos_omega = 10;	// <= _N
	
	// ----------- declaration of matrices ---------------
	Matrix<double, _n , _n> A_arg;		// n x n
	Matrix<double, _n, _m> B_arg;		// n x m; resizng for non-square matrices doesn't work
	Matrix<double, _n, 1> s_arg;		// n x 1
	Matrix<double, _n, _n> Q_tilde_arg;	// n x n
	Matrix<double, _n, _n> Q_tilde_f_arg;	// n x n
	Matrix<double, _m, _m> R_arg;			// m x m
	Matrix<double, _m, _n> K_arg;			// m x n
	Matrix<double, _nSt, _n> Fx_arg[_N];		// _nSt x n, [_N]
	Matrix<double, _nSt, 1> fx_arg[_N];		// _nSt x 1, [_N] 
	Matrix<double, _nInp, _m> Fu_arg[_N];		// _nInp x m, [_N]
	Matrix<double, _nInp, 1> fu_arg[_N];		// _nInp x 1, [_N]
	
	Matrix<double, _nF_xTheta, _n> F_xTheta_arg;	// _nF_xTheta x n
	Matrix<double, _nF_xTheta, _m> F_theta_arg;	// _nF_xTheta x m
	Matrix<double, _nF_xTheta, 1> f_xTheta_arg;	// _nF_xTheta x 1
	
	Matrix<double, _n, _n> Lm_arg;		// n x n
	Matrix<double, _n, _m> Mm_arg;		// n x m
	Matrix<double, _n, 1> tm_arg;		// n x 1
	Matrix<double, _n, 1> x_hat_arg;		// n x 1, state estimate
	Matrix<double, _n, 1> x_star_arg[_N];	// n x 1, [_N], tracking
	Matrix<double, _m, 1> u_opt;			// m x 1, optimal input is saved there
	
	
	// ---------- no changes necessary ---------
	double kappa_arg;	// for PhaseII
	double kappa_PhaseI_arg;	// for PhaseI	
	int n_iter_arg;
	double mu_arg;
	double eps_barrier_arg;
	double eps_nt_arg;
	double eps_normRp_arg;
	double eps_ls_arg;	
	double alpha_ls_arg;
	double beta_ls_arg;
	double reg_arg;		// regularization term for PhaseII
	double reg_PhaseI_arg;	// regularization term for PhaseI
	
	// ---------------- read from binary file ------------------
	ifstream fin;			 	// Definition input file object
	fin.open("ConstrParam.bin", ios::binary);	//	open file
	if (!fin.is_open())
	{
		cout << "File open error \n"; 
		return 1;
	}
	
	// read 
	fin.read((char *) &kappa_arg, sizeof(double));
	fin.read((char *) &kappa_PhaseI_arg, sizeof(double));
	fin.read((char *) &n_iter_arg, sizeof(int));
	fin.read((char *) &mu_arg, sizeof(double));
	fin.read((char *) &eps_barrier_arg, sizeof(double));
	fin.read((char *) &eps_nt_arg, sizeof(double));
	fin.read((char *) &eps_normRp_arg, sizeof(double));
	fin.read((char *) &eps_ls_arg, sizeof(double));
	fin.read((char *) &alpha_ls_arg, sizeof(double));
	fin.read((char *) &beta_ls_arg, sizeof(double));
	fin.read((char *) &reg_arg, sizeof(double));
	fin.read((char *) &reg_PhaseI_arg, sizeof(double));
	

	// read A_arg
	for (int i = 0; i <= _n-1; i++)		
	{
		for (int j = 0; j <= _n-1; j++)
		{
			fin.read((char *) &A_arg(j,i), sizeof(double));
		}
	}

	// read B_arg
	for (int i = 0; i <= _m-1; i++)	// #columns
	{
		for (int j = 0; j <= _n-1; j++) // #rows
		{
			fin.read((char *) &B_arg(j,i), sizeof(double));
		}
	}
	
	// read s_arg
	for (int i = 0; i <= _n-1; i++)	// #columns
	{	
		fin.read((char *) &s_arg(i,0), sizeof(double));
	}

	// read Q_tilde_arg
	for (int i = 0; i <= _n-1; i++)	// #columns
	{
		for (int j = 0; j <= _n-1; j++) // #rows
		{
			fin.read((char *) &Q_tilde_arg(j,i), sizeof(double));
		}
	}
	
	// read Q_tilde_f_arg
	for (int i = 0; i <= _n-1; i++)	// #columns
	{
		for (int j = 0; j <= _n-1; j++) // #rows
		{
			fin.read((char *) &Q_tilde_f_arg(j,i), sizeof(double));
		}
	}
	
	// read R_arg
	for (int i = 0; i <= _m-1; i++)	// #columns
	{
		for (int j = 0; j <= _m-1; j++) // #rows
		{
			fin.read((char *) &R_arg(j,i), sizeof(double));
		}
	}
	
	// read Fx_arg[]
	for (int k = 0; k <= _N-1; k++)
	{
		for (int i = 0; i <= _n-1; i++)	// #columns
		{
			for (int j = 0; j <= _nSt-1; j++) // #rows
			{
				fin.read((char *) &Fx_arg[k](j,i), sizeof(double));
			}
		}
	}
	
	// read fx_arg[]
	for (int k = 0; k <= _N-1; k++)
	{
		for (int i = 0; i <= _nSt-1; i++)	// #columns
		{
				fin.read((char *) &fx_arg[k](i,0), sizeof(double));
		}
	}
	
	// read Fu_arg[]
	for (int k = 0; k <= _N-1; k++)
	{
		for (int i = 0; i <= _m-1; i++)	// #columns
		{
			for (int j = 0; j <= _nInp-1; j++) // #rows
			{
				fin.read((char *) &Fu_arg[k](j,i), sizeof(double));
			}
		}
	}
	
	// read fu_arg[]
	for (int k = 0; k <= _N-1; k++)
	{
		for (int i = 0; i <= _nInp-1; i++)	// #columns
		{
				fin.read((char *) &fu_arg[k](i,0), sizeof(double));
		}
	}
	
	// read F_xTheta_arg
	for (int i = 0; i <= _n-1; i++)	// #columns
	{
		for (int j = 0; j <= _nF_xTheta-1; j++) // #rows
		{
			fin.read((char *) &F_xTheta_arg(j,i), sizeof(double));
		}
	}
	
	// read F_theta_arg
	for (int i = 0; i <= _m-1; i++)	// #columns
	{
		for (int j = 0; j <= _nF_xTheta-1; j++) // #rows
		{
			fin.read((char *) &F_theta_arg(j,i), sizeof(double));
		}
	}
	
	// read f_xTheta_arg
	for (int i = 0; i <= _nF_xTheta-1; i++)	// #columns
	{	
		fin.read((char *) &f_xTheta_arg(i,0), sizeof(double));
	}

	// read K_arg
	for (int i = 0; i <= _n-1; i++)	// #columns
	{
		for (int j = 0; j <= _m-1; j++) // #rows
		{
			fin.read((char *) &K_arg(j,i), sizeof(double));
		}
	}

	fin.close();				// close file
	
	/*
	cout << "kappa_arg: " << kappa_arg << endl << endl;
	cout << "kappa_PhaseI_arg " << kappa_PhaseI_arg << endl << endl;
	cout << "n_iter_arg: " << n_iter_arg << endl << endl;
	cout << "mu_arg: " << mu_arg << endl << endl;
	cout << "eps_barrier_arg: " << eps_barrier_arg  << endl << endl;
	cout << "eps_nt_arg: " << eps_nt_arg << endl << endl;
	cout << "eps_normRp_arg: " << eps_normRp_arg << endl << endl;
	cout << "eps_ls_arg: " << eps_ls_arg << endl << endl;
	cout << "alpha_ls_arg: " << alpha_ls_arg << endl << endl;
	cout << "beta_ls_arg: " << beta_ls_arg << endl << endl;
	cout << "n_iter_arg: " << n_iter_arg << endl << endl;
	cout << "reg_arg: " << reg_arg << endl << endl;
	cout << "reg_PhaseI_arg: " << reg_PhaseI_arg << endl << endl;
	cout << "A_arg: " << endl << A_arg << endl << endl;
	cout << "B_arg: " << endl << B_arg << endl << endl;
	cout << "s_arg:" << endl << s_arg << endl << endl;
	cout << "Q_tilde_arg: " << endl << Q_tilde_arg << endl << endl;
	cout << "Q_tilde_f_arg: " << endl << Q_tilde_f_arg << endl << endl;
	cout << "R_arg: " << endl << R_arg << endl << endl;
	
	for (int i = 0; i<= _N-1; i++)
	{
		cout << "Fx_arg[" << i << "]: " << endl << Fx_arg[i] << endl << endl;
	}
	for (int i = 0; i <= _N-1; i++)
	{
		cout << "fx_arg[" << i << "]: " << endl << fx_arg[i] << endl << endl;
	}
	for (int i = 0; i <= _N-1; i++)
	{
		cout << "Fu_arg[" << i << "]: " << endl << Fu_arg[i] << endl << endl;
	}
	for (int i = 0; i <= _N-1; i++)
	{
		cout << "fu_arg[" << i << "]: " << endl << fu_arg[i] << endl << endl;
	}
	cout << "F_xTheta_arg: " << endl << F_xTheta_arg << endl << endl;
	cout << "F_theta_arg: " << endl << F_theta_arg << endl << endl;
	cout << "f_xTheta_arg: " << endl << f_xTheta_arg << endl << endl;
	cout << "K_arg:" << endl << K_arg << endl << endl;
	*/
	
	
	// -------- object instantiation -----------
	 LBmpcTP<double, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> myObj(				// constructor
	 		 	  	 			 	kappa_arg, kappa_PhaseI_arg, n_iter_arg,  mu_arg,  eps_barrier_arg,  eps_nt_arg, eps_normRp_arg, eps_ls_arg, 
									alpha_ls_arg, beta_ls_arg, reg_arg, reg_PhaseI_arg, A_arg, B_arg, Q_tilde_arg, Q_tilde_f_arg, R_arg, Fx_arg, 
	 		 	 	 	 		   	fx_arg, Fu_arg, fu_arg, F_xTheta_arg, F_theta_arg, f_xTheta_arg, K_arg, s_arg);
		 	


	
	// ----------- SPECIFY arguments for step() ---------------
	// -------- they are updated at each time step ------------
	
	Lm_arg << 1,  2,  0,  1,  2,
	      -2,  1.3,  1,  2,  2.3,
	       1,  2, -1,  0,  -1,
	       1,  2,  2,  -2.3,  1,
	      0, 0,  2,  1.4, -2;
	
	Mm_arg << 1,  1.4,
	   2,  -1,
	   1,   2,
	   0,   0,
	   2, -1;
	
	tm_arg << -1, 2, 1, 1, 2;
	x_hat_arg << 3, 3, -2, 3, 4;

	srand((unsigned)time(0));
	
	for(int i = 0; i <= _N-1; i++)
	{
		x_star_arg[i].setZero();
	}
	
	// the following matrices are only needed for test purposes
	Matrix<double, _N*(_m+_n+_n)+_m,1> z_warm_arg;	// N*(m+n+n)+m	
	
	z_warm_arg.setRandom();
	z_warm_arg = 10*z_warm_arg;


	u_opt = myObj.step(	Lm_arg, Mm_arg, tm_arg,	x_hat_arg, z_warm_arg, x_star_arg);
	cout << "optimal input:" << endl << u_opt << endl << endl;

	return 0;
	
}



















