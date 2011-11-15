// mainLBmpcTP.cpp
// example file to test simple examples 
// date: November 08, 2011
// author: Xiaojing ZHANG

// matrices are imported from binary file created by MATLAB


#include <iostream>		// I-O
#include <string>
#include <Eigen/Dense>	// matrix computation
#include "LBmpcTP.h"	// class template
#include <time.h>

using namespace Eigen;
using namespace std;

int main(int argc, const char* argv[])
{
  // ------ SPECIFY parameters -------
  const int _N = 60; // MPC horizon
  const int _m = 3; // #input
  const int _n = 10; // #states
  const int _nSt = 20; // # state constraints
  const int _nInp = 6; // # state constraints
  const int _nF_xTheta = 72; // # Omega constraints
  const int _pos_omega = 1; // <= _N
  const string default_fileName = "ConstrParam.bin";
  string filename;
  bool verbose = 0; // '0' if it should shut up
  int status; // stores status code

  if (argc == 1)
  {
    filename = default_fileName;
  }
  else if (argc == 2)
  {
    filename = argv[1];
  }
  else
  {
    cerr << "Usage: " << argv[0] << " [path-to-ConstrParam.bin]" << endl;
    return 1;
  }

#ifdef EIGEN_VECTORIZE
  cout << "EIGEN_VECTORIZE set!" << endl;
#else
  cout << "EIGEN_VECTORIZE not set!" << endl;
#endif
#ifdef NDEBUG
  cout << "NDEBUG set" << endl;
#endif

  // -------- object instantiation -----------
  // fileName contains name of file with parameters
  // bool verbose: 1 or 0
  LBmpcTP<double, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> myObj(filename.c_str(), verbose);

  Matrix<double, _n, _n> A;
  Matrix<double, _n, _m> B;
  Matrix<double, _n, 1> s;
  myObj.get_A(A);
  myObj.get_B(B);
  myObj.get_s(s);

  // ----------- SPECIFY arguments for step() ---------------
  // those matrices are computed externely
  // -----------  sizes of matrices ---------------
  Matrix<double, _n, _n> Lm; // n x n
  Matrix<double, _n, _m> Mm; // n x m
  Matrix<double, _n, 1> tm; // n x 1
  Matrix<double, _n, 1> x_hat; // n x 1, state estimate
  Matrix<double, _n, 1> x_star[_N]; // n x 1, [_N], tracking
  Matrix<double, _m, 1> u_opt; // m x 1, optimal input is saved there

  // -------- they are updated at each time step ------------
  Lm.setZero();
  //Lm << 1, 2, 0, 1, 2, -2, 1.3, 1, 2, 2.3, 1, 2, -1, 0, -1, 1, 2, 2, -2.3, 1, 0, 0, 2, 1.4, -2;
  Mm.setZero();
//  Mm << 1, 1.4, 2, -1, 1, 2, 0, 0, 2, -1;

  tm.setZero();
//  tm << -1, 2, 1, 1, 2;
//  x_hat << 3, 3, -2, 3, 4;
  x_hat << 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0;

  for (int i = 0; i <= _N - 1; i++)
  {
    x_star[i].setZero();
    x_star[i][8] = -2.0;
  }
  // x_star_arg[0] << 29.1867 ,  20.3181,   36.6838,   10.5584,  -24.6923;
  // x_star_arg[1] << 401.2845,  262.2318,  -73.8211, -285.3312,  391.4877;
  // x_star_arg[2] << -1.4669*1000  ,  0.0849*1000 ,   0.1898*1000 ,   2.8506*1000 ,  -1.8977*1000;
  // x_star_arg[3] << -0.8542*1000  ,  9.2976*1000   , 7.0920*1000 ,  -1.5346*1000 ,   4.6926*1000;


  for (int j = 0; j < 1500; j++)
  {
    timespec start_rt;
    timespec end_rt;
    clock_gettime(CLOCK_REALTIME, &start_rt);
    status = myObj.step(Lm, Mm, tm, x_hat, x_star);
    clock_gettime(CLOCK_REALTIME, &end_rt);
    double duration = (end_rt.tv_sec - start_rt.tv_sec) + 1e-9*(end_rt.tv_nsec - start_rt.tv_nsec);
    cout << "status: " << status << " dur: " << duration;
    cout << " iter: " << myObj.n_iter_last;
    if (!status)
    {
      u_opt = myObj.u_opt;
      cout << " u_opt: " << u_opt.transpose();
      cout << " x_hat: " << x_hat.transpose();
      cout << endl;
    } else
    {
      cout << endl;
    }

    x_hat = A * x_hat + B * u_opt + s;
  }


  return 0;
}

