// mainLBmpcTP.cpp
// example file to test simple examples 
// date: November 08, 2011
// author: Xiaojing ZHANG

// matrices are imported from binary file created by MATLAB


#include <iostream>		// I-O
#include <string>
#include <Eigen/Dense>	// matrix computation
#include "LBmpcTP.h"	// class template
#include <time.h>	// for Linux
// #include <sys/time.h> 
#include <sys/time.h>

using namespace Eigen;
using namespace std;

const int steps = 1000; // number of simulations steps

template<int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  class LBMPCAxis
  {
  private:
    LBmpcTP<double, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>* myObj;
    Matrix<double, _n, _n> A;
    Matrix<double, _n, _m> B;
    Matrix<double, _n, 1> s;
    Matrix<double, _n, _n> Lm; // n x n
    Matrix<double, _n, _m> Mm; // n x m
    Matrix<double, _n, 1> tm; // n x 1
    Matrix<double, _n, 1> x_hat; // n x 1, state estimate
    Matrix<double, _n, 1> x_star[_N]; // n x 1, [_N], tracking
    Matrix<double, _m, 1> u_opt; // m x 1, optimal input is saved there
    // ------ SPECIFY parameters -------
    // only works for _N <= 243, other George cannot execute main file:
    //  const int _N = 15; // MPC horizon
    //  const int _m = 3; // #input
    //  const int _n = 10; // #states
    //  const int _nSt = 20; // # state constraints
    //  const int _nInp = 6; // # input constraints
    //  const int _nF_xTheta = 36; // # Omega constraints
    //  const int _pos_omega = 1; // <= _N
    //  const string default_fileName = "ConstrParam.bin";
    string filename;
    bool verbose; // '0' if it should shut up
    int status; // stores status code

    // add some noise
    // srand ( 23 );
    // Matrix<double, _n, steps> noise;
    // noise.setRandom();
    // cout << "noise: " << noise << endl;
    // noise = 0*0.05*noise;
  public:
    LBMPCAxis(string filename)
    {

      cout << endl << endl;
#ifdef EIGEN_VECTORIZE
      cout << "EIGEN_VECTORIZE set!" << endl;
#else
      cout << "EIGEN_VECTORIZE not set!" << endl;
#endif
#ifdef NDEBUG
      cout << "NDEBUG set" << endl;
#endif
      verbose = 0;
      myObj = new LBmpcTP<double, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> (filename.c_str(), verbose);

      // -------- object instantiation -----------
      // fileName contains name of file with parameters
      // bool verbose: 1 or 0
      myObj->get_A(A);
      myObj->get_B(B);
      myObj->get_s(s);

      // ----------- SPECIFY arguments for step() ---------------
      // those matrices are computed externely
      // -----------  sizes of matrices ---------------

      // -------- they are updated at each time step ------------
      Lm.setZero();
      Mm.setZero();

      tm.setZero();

      initialize();
    }

    void initialize()
    {
      x_hat.setZero();
      //x_hat << 0, 0, 0, 0, 0, 0, 0, 0, -2, 0;
      // x_hat << 0,0,0,0,0,0,0,0,-1,0;
      double zcmd = -2.0;

    }

    void step()
    {
      status = myObj->step(Lm, Mm, tm, x_hat, x_star);
      if (!status)
      {
        u_opt = myObj->u_opt;
        x_hat = A * x_hat + B * u_opt + s;
      }
      else
      {
        cerr << "error code: " << status << endl;
        cerr << "starting position was: " << x_hat.transpose() << endl;
        cout << "most optimal u_opt: " << myObj->u_opt.transpose() << endl;

      }
    }

    void set_x_hat(const Matrix<double, _n, 1>& x_hat_set)
    {
      x_hat = x_hat_set;
    }

    void set_x_star(const Matrix<double, _n, 1>& x_star_set)
    {
      for (int i = 0; i <= _N - 1; i++)
      {
        x_star[i] = x_star_set;
      }
    }

    void get_x_hat(Matrix<double, _n, 1>& x_hat_out)
    {
      x_hat_out = x_hat;
    }

    void get_u_opt(Matrix<double, _m, 1>& u_opt_out)
    {
      u_opt_out = u_opt;
    }

    int get_niter_last()
    {
      return myObj->n_iter_last;
    }

  }; // class

int main(int argc, const char* argv[])
{

  const int _N = N_MPC_STEPS;
  const int nxy = 4;
  const int m = 1;
  const int nz = 2;
  int numIter = 0;
  const int nF_xTheta_xy = 14; // need to manually input this from Matlab output
  const int nF_xTheta_z = 8; // need to manually input this from Matlab output

  struct timeval start; // added by George
  struct timeval end; // ditto
  double elapsedTime = 0; // ditto
  double timeTmp; // ditto
  LBMPCAxis<nxy, m, _N, nxy * 2, m * 2, nF_xTheta_xy, 1> lbmpc_x(argv[1]);
  LBMPCAxis<nxy, m, _N, nxy * 2, m * 2, nF_xTheta_xy, 1> lbmpc_y(argv[1]);
  LBMPCAxis<nz, m, _N, nz * 2, m * 2, nF_xTheta_z, 1> lbmpc_z(argv[2]);

  Matrix<double, nxy, 1> x_hatx, x_starx;
  Matrix<double, nxy, 1> x_haty, x_stary;
  Matrix<double, nz, 1> x_hatz, x_starz;

  // Initial condition:
  x_hatx << 0, 0, 0, 0;
  x_haty << 0, 0, 0, 0;
  x_hatz << -2, 0;

  // Desired state:
  x_starx << 0, 0, 0, 0;
  x_stary << 0, 0, 0, 0;
  x_starz << -1, 0;

  lbmpc_x.set_x_hat(x_hatx);
  lbmpc_y.set_x_hat(x_haty);
  lbmpc_z.set_x_hat(x_hatz);
  lbmpc_x.set_x_star(x_starx);
  lbmpc_y.set_x_star(x_stary);
  lbmpc_z.set_x_star(x_starz);

  //Matrix<double, nxy * 2, 1> x_hat;

  Matrix<double, m, 1> u_optx, u_opty, u_optz;

  lbmpc_x.get_x_hat(x_hatx);
  lbmpc_y.get_x_hat(x_haty);
  lbmpc_z.get_x_hat(x_hatz);

  cout << "x_hatx = " << x_hatx.transpose() << endl;
  cout << "x_haty = " << x_haty.transpose() << endl;

  timespec start_rt; // linux
  timespec end_rt; // linux
  double duration;

  for (int j = 0; j < steps; j++)
  {
    //clock_gettime(CLOCK_REALTIME, &start_rt); // linux
    gettimeofday(&start, NULL); // George's Mac

    lbmpc_x.step();
    lbmpc_y.step();
    lbmpc_z.step();
    //clock_gettime(CLOCK_REALTIME, &end_rt); // linux
    gettimeofday(&end, NULL); // George's Mac
    //duration = (end_rt.tv_sec - start_rt.tv_sec) + 1e-9 * (end_rt.tv_nsec - start_rt.tv_nsec); // linux
    duration = (end.tv_sec * 1000 + end.tv_usec / 1000) - (start.tv_sec * 1000 + start.tv_usec / 1000);
    elapsedTime += duration;


    lbmpc_x.get_x_hat(x_hatx);
    lbmpc_y.get_x_hat(x_haty);
    lbmpc_z.get_x_hat(x_hatz);
    lbmpc_x.get_u_opt(u_optx);
    lbmpc_y.get_u_opt(u_opty);
    lbmpc_z.get_u_opt(u_optz);
    cout << "j = " << j << endl;
    cout << "X: x = " << x_hatx.transpose() << " u_opt = " << u_optx << " iter = " << lbmpc_x.get_niter_last() << endl;
    cout << "Y: x = " << x_haty.transpose() << " u_opt = " << u_opty << " iter = " << lbmpc_y.get_niter_last() << endl;
    cout << "Z: x = " << x_hatz.transpose() << " u_opt = " << u_optz << " iter = " << lbmpc_z.get_niter_last() << endl;
    cout << "dur = " << duration << " [ms]" << endl;
    // x_hat = A * x_hat + B * u_opt + s + noise.col(j);
  }
  cout << endl;
  cout << "Finished " << steps << " steps in " << elapsedTime << " ms" << endl;
  cout << (steps / elapsedTime) << " steps/ms, " << (elapsedTime / steps) << " ms/step" << endl;

  return 0;
}

