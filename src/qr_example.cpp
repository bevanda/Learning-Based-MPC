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
  const int _N = N_MPC_STEPS; // MPC horizon
  const int _m = 3; // #input
  const int _n = 10; // #states
  const int _nSt = 20; // # state constraints
  const int _nInp = 6; // # input constraints
  const int _nF_xTheta = 36; // # Omega constraints
  const int _pos_omega = 1; // <= _N
  const string default_fileName = "ConstrParam.bin";
  string filename;
  bool verbose = 0; // '0' if it should shut up
  int status; // stores status code

  // add some noise
  // srand ( 23 );
  const int steps = 1000; // number of simulations steps
  // Matrix<double, _n, steps> noise;
  // noise.setRandom();
  // cout << "noise: " << noise << endl; 
  // noise = 0*0.05*noise;

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

  // Use some small constant values to have a nonzero oracle update...
  Lm.setZero();
  Lm(3, 2) = Lm(3, 3) = Lm(7, 6) = Lm(7, 7) = 0.01;
  Mm.setZero();
  Mm(3,0) = Mm(7,1) = 0.001;
  Mm(9,2) = 1e-5;
  tm.setZero();
  tm[1] = tm[3] = tm[5] = tm[7] = 1e-3;
  tm[8] = 1e-5;
  tm[9] = Mm(9,2) * 80;

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
    //    if (j == 750)
    //    {
    //      for (int i = 0; i <= _N - 1; i++)
    //      {
    //        x_star[i].setZero();
    //        x_star[i][0] = 1.0;
    //        x_star[i][8] = -3.0;
    //      }
    //    }
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
  cout << (steps / elapsedTime) << " steps/ms, " << (elapsedTime / steps) << " ms/step" << endl;

  return 0;
}

