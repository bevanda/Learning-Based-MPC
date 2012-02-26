// mainLBmpcTP.cpp
// example file to test simple examples 
// date: November 08, 2011
// author: Xiaojing ZHANG

// matrices are imported from binary file created by MATLAB


#include <iostream>		// I-O
#include <string>
#include <Eigen/Dense>	// matrix computation
#include "LBmpcTP.h"	// class template
// #include <time.h>	// for Linux
// #include <sys/time.h> 
#include <sys/time.h>

using namespace Eigen;
using namespace std;

int main(int argc, const char* argv[])
{
  // ------ SPECIFY parameters -------
  // only works for _N <= 243, other George cannot execute main file:
  const int _N = 15; // MPC horizon
  //const int _m[3] = {1,1,1}; // #input
  const int _n = 10; // #states
  const int _nSt = 20; // # state constraints
  const int _nInp = 6; // # input constraints
  const int _nF_xTheta = 36; // # Omega constraints
  //const string default_fileName = "ConstrParam.bin";
  string filename_x, filename_y, filename_z;
  bool verbose = 0; // '0' if it should shut up
  int status; // stores status code

  // add some noise
  // srand ( 23 );
  const int steps = 1000; // number of simulations steps
  // Matrix<double, _n, steps> noise;
  // noise.setRandom();
  // cout << "noise: " << noise << endl; 
  // noise = 0*0.05*noise;

  filename_x = argv[1];
  filename_y = argv[2];
  filename_z = argv[3];

  cout << endl << endl;
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
  LBmpcTP<double, 4, 1, _N, 8, 2, 0, 1> lbmpctp_x(filename_x.c_str(), verbose);
  LBmpcTP<double, 4, 1, _N, 8, 2, 0, 1> lbmpctp_y(filename_y.c_str(), verbose);
  LBmpcTP<double, 2, 1, _N, 4, 2, 0, 1> lbmpctp_y(filename_y.c_str(), verbose);
  Matrix<double, 4, 4> Ax, Ay;
  Matrix<double, 2, 2> Az, Az;
  Matrix<double, 4, 1> Bx, By;
  Matrix<double, 4, 1> Bz;
  Matrix<double, 4, 1> sx, sy;
  Matrix<double, 2, 1> sz;
  lbmpctp_x.get_A(Ax);
  lbmpctp_y.get_A(Ay);
  lbmpctp_z.get_A(Az);
  lbmpctp_x.get_B(Bx);
  lbmpctp_y.get_B(By);
  lbmpctp_z.get_B(Bz);
  lbmpctp_x.get_s(sx);
  lbmpctp_y.get_s(sy);
  lbmpctp_z.get_s(sz);

  // ----------- SPECIFY arguments for step() ---------------
  // those matrices are computed externely
  // -----------  sizes of matrices ---------------
  Matrix<double, 4, 4> Lmx, Lmy; // n x n
  Matrix<double, 2, 2> Lmz; // n x n
  Matrix<double, 4, 1> Mmx, Mmy; // n x m
  Matrix<double, 2, 1> Mmz; // n x m
  Matrix<double, 4, 1> tmx, tmy; // n x 1
  Matrix<double, 2, 1> tmz; // n x 1
  Matrix<double, 4, 1> x_hatx, x_haty; // n x 1, state estimate
  Matrix<double, 2, 1> x_hatz, x_hatz; // n x 1, state estimate
  Matrix<double, _n, 1> x_star[_N]; // n x 1, [_N], tracking
  Matrix<double, _m[0], 1> u_opt; // m x 1, optimal input is saved there

  // -------- they are updated at each time step ------------
  Lm.setZero();
  Mm.setZero();

  tm.setZero();
  x_hat << 0, 0, 0, 0, 0, 0, 0, 0, -2, 0;
	// x_hat << 0,0,0,0,0,0,0,0,-1,0;
  double zcmd = -2.0;

  for (int i = 0; i <= _N - 1; i++)
  {
    x_star[i].setZero();
    x_star[i][8] = -1.0;
	// x_star[i][0] = 1.0;
    // x_star[i][8] = -3.0;
  }
  int numIter = 0;

  struct timeval start; // added by George
  struct timeval end; // ditto
  double elapsedTime = 0; // ditto
  double timeTmp; // ditto

  for (int j = 0; j < steps; j++)
  {

    // change alt. setpoint after awhile:
    if (j == 750)
    {
      for (int i = 0; i <= _N - 1; i++)
      {
        x_star[i].setZero();
        x_star[i][0] = 1.0;
        x_star[i][8] = -3.0;
      }
    }
    // timespec start_rt;	// linux
    // timespec end_rt;		// linux
    // clock_gettime(CLOCK_REALTIME, &start_rt);	// linux
    gettimeofday(&start, NULL); // George's Mac
    status = myObj.step(Lm, Mm, tm, x_hat, x_star);
    gettimeofday(&end, NULL); // George's Mac
    timeTmp = (end.tv_sec * 1000 + end.tv_usec / 1000) - (start.tv_sec * 1000 + start.tv_usec / 1000); // George's Mac
    elapsedTime += timeTmp; // George's Mac

    // clock_gettime(CLOCK_REALTIME, &end_rt);		// linux
    // double duration = (end_rt.tv_sec - start_rt.tv_sec) + 1e-9*(end_rt.tv_nsec - start_rt.tv_nsec);	// linux
    // cout << "step: " << j ;
    // cout << " status: " << status;
    // cout << " iter: " << myObj.n_iter_last;
    numIter = numIter + myObj.n_iter_last;
 
    if (!status)
    {
		cout << setprecision(8);
      u_opt = myObj.u_opt;
      cout << j;
      cout << " iter: " << myObj.n_iter_last;
      cout << " u_opt: " << u_opt.transpose();
      cout << " | x_hat: " << x_hat.transpose();
      cout << endl << endl;
    }
    else
    {
      cerr << "error code: " << status << " at step: " << j << endl;
      cerr << "starting position was: " << x_hat.transpose() << endl;
      cout << "most optimal u_opt: " << myObj.u_opt.transpose() << endl;
      // return 1;
      u_opt.setZero();
	return 1;
    }

    // x_hat = A * x_hat + B * u_opt + s + noise.col(j);
    x_hat = A * x_hat + B * u_opt + s;
  }

  cout << endl << endl;
  cout << "STAT: " << endl;
  cout << "number of Newton iterations: " << numIter << endl;
  cout << "time elapsed: " << elapsedTime << " [ms]" << endl;
  cout << "average: " << elapsedTime / numIter << " [ms/step]" << endl;
  return 0;
}

