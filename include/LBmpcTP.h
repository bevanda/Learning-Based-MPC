// Ompctp.h -- a Learning Based MPC class template
// header for class LBmpcTP
// date: November 8, 2011
// author: Xiaojing ZHANG
// version: 0.2


/* Remarks:

 2) horizon >= 3 assumed
 3) issues with LLT, L_diag[2*Horizon] for both Y and Phi --> save it into normal matrices
 4) issues of defining length and Type for vector d, use diagonal Matrix for d^2
 5) Optimize index calculations, merge different loops into 1 loop
 6) avoid using namespace Eigen/std
 11) compGamma(): use vector addition rather than piece by piece addition?
 12) directly use _arg from arguments rather than copying, esp. during step();


 Improvements:
 1) different step length for primal and dual step (Numerical Optimization Nocedal/Wright)
 2) different/better termination condition?
 3) adaptive damping factor
 */

#ifndef LBMPCTP_H_
#define LBMPCTP_H_

#include <iostream>
#include <fstream>		// read binary data
#include <Eigen/Dense>
#include <iomanip>		// for setprecision (9)
#include <cmath>

using namespace Eigen;
using namespace std;

/* --------- input parameters for template -----------
 Type:		double or float as matrix elements
 _n:			dimension of state
 _m:			dimension input
 _N:			dimension MPC horizon
 _nSt:		number of state constraints (>= n for full-rank)
 _nInp:		number of input constraints (>= m for full-rank)
 _nF_xTheta:	number of constraints of F_xTheta*x[m+N]+F_theta*theta <= f_xTheta
 _pos_omega: (x_bar[.+_pos_omega] , theta) \in Omega
 */
const int Nhoriz = 60;

template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  class LBmpcTP
  {
  private:
    // ----------- input variables ------------
    int n_iter; // number of Newton iterations for fixed kappa in PhaseII
    double alpha; // step size
    int offset; // = _n + _n + _m: often used variable
    double resNorm_H;
    double resNorm_C;
    double resNorm_P;
    double muNorm;
    double damp;
    double reg; // regularization term
    double compTime; // stores computational time
    bool verbose; // if set, program will print status updates

    Matrix<Type, _n, _n> A; // dynamics matrix A, n x n
    Matrix<Type, _n, _n> A_transp; // transpose of above
    Matrix<Type, _n, _m> B; // input matrix B, n x m
    Matrix<Type, _m, _n> B_transp; // transpose of above
    Matrix<Type, _n, 1> s; // affine offsets in linear dynamics
    Matrix<Type, _n, _n> Q_tilde; // cost matrix on state x_tilde, p.d., n x n
    Matrix<Type, _n, _n> Q_tilde_f; // final cost matrix, p.d.
    Matrix<Type, _m, _m> R; // cost matrix on input, p.d. if Fu not full rank
    Matrix<Type, _n, 1> q_tilde_vec[_N]; // cost matrix
    Matrix<Type, _m, 1> r_vec[_N];
    Matrix<Type, _nSt, _n> Fx[_N]; // array of full-rank state constraint matrices
    Matrix<Type, _n, _nSt> Fx_transp[_N]; // transpose of above
    Matrix<Type, _nSt, 1> fx[_N]; // array containing right-hand-side of state constraints
    Matrix<Type, _nInp, _m> Fu[_N]; // array of (full-rank) input constraint matrices
    Matrix<Type, _m, _nInp> Fu_transp[_N]; //transpose of above
    Matrix<Type, _nInp, 1> fu[_N]; // array containing right-hand-side of input constraints
    Matrix<Type, _nF_xTheta, _n> F_xTheta; // F_xTheta*x[m+j]+F_theta*theta <= f_xTheta
    Matrix<Type, _n, _nF_xTheta> F_xTheta_transp; // transpose of above
    Matrix<Type, _nF_xTheta, _m> F_theta; // full-rank constraint matrix on theta
    Matrix<Type, _m, _nF_xTheta> F_theta_transp; // transpose of above
    Matrix<Type, _nF_xTheta, 1> f_xTheta; // right-hand-side of constraint above

    Matrix<Type, _m, _n> K; // Gain matrix, s.t. (A+B*K) is stable/Schur
    Matrix<Type, _n, _m> K_transp;

    const Matrix<Type, _n, 1> *x_hat; // state estimate of current state x_hat[m]

    // search directions dz and dnu for all variables z, nu, lambda, t
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> z; // z(+) = z + alpha*dz
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> dz;
    Matrix<Type, 2 * _N * _n, 1> nu; // nu(+) = nu + alpha*dnu
    Matrix<Type, 2 * _N * _n, 1> dnu;
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> lambda; // > 0
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> dlambda;
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> slack; // > 0
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> dslack;

    // residua
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> r_H;
    Matrix<Type, 2 * _N * _n, 1> r_C;
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> r_P;
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> r_T;
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> r_d_bar;

    // ----------- auxiliary variables ------------
    Matrix<Type, _n, 1> tm_tilde; // tm_tilde = s + tm
    Matrix<Type, _n, _m> S; // S = K'*R
    Matrix<Type, _m, _n> S_transp; // S_trans = S'
    Matrix<Type, _n, _n> Q_bar; // Q_bar = K'*R*K = S*K
    Matrix<Type, _n, 1> q_bar_vec[_N];
    Matrix<Type, _n, _n> A_bar; // A_bar = A+B*K
    Matrix<Type, _n, _n> A_bar_transp; // transpose of above
    Matrix<Type, _n, _n> Am_tilde; // Am_tilde = A + Lm
    Matrix<Type, _n, _n> Am_tilde_transp; // transpose of above
    Matrix<Type, _n, _m> Bm; // Bm = B + Mm
    Matrix<Type, _m, _n> Bm_transp; // transpose of above
    Matrix<Type, _n, _n> Bm_bar; // Bm_bar = Bm * K
    Matrix<Type, _n, _n> Bm_bar_transp; // transpose of above
    Matrix<Type, _nInp, _n> Fu_bar[_N]; // array of (full-rank) input constraint matrices, #row may change, but <= nInp
    Matrix<Type, _n, _nInp> Fu_bar_transp[_N]; //transpose of above

    // used in compResiduals(),
    Matrix<Type, _m, 1> c_tmp; // use this to temporarily copy c-blocks out of z
    Matrix<Type, _n, 1> x_bar_tmp; // use this to temporarily copy x_bar-blocks out of z

    // Variables needed to represent Phi -> see documentation
    Matrix<Type, _m, _m> Omicron[_N + 1];
    Matrix<Type, _n, _n> Rho[_N];
    Matrix<Type, _n, _m> Sigma[_N];

    // Variable to represent Y = C*Phi*C'
    Matrix<Type, _n, _n> Y[3][2 * _N]; // Y[i][j] = Y_{i+1,j+1}, i=0,1,2
    // Y[1][0], Y[2][0], Y[2][1] not used

    // !!!!!!! Variables to represent L: L*L'=Y
    // !!!!!!! need explicitly: double, Dynamic, array-length, e.g. [2*Horizon]
    LLT<Matrix<double, Dynamic, Dynamic> , Lower> L_diag[2 * Nhoriz]; // ****** diagonal matrices of L computed using Cholesky
    Matrix<Type, _n, _n> L_diag_transp[2 * _N];
    Matrix<Type, _n, _n> L_offDiag[2][2 * _N]; // off-diag matrices are general square matrices
    Matrix<Type, _n, _n> L_offDiag_transp[2][2 * _N]; // transpose of above

    // !!!!!!!! Variables to represent Phi = L_Phi*L_Phi'
    // !!!!!!!! need explicitly: double, Dynamic, array-length, e.g. horizon = 50
    LLT<Matrix<double, Dynamic, Dynamic> , Lower> LOmicron_diag[Nhoriz + 1]; //*** [horizon +1]
    Matrix<double, _m, _m> LOmicron_diag_transp[_N + 1];
    LLT<Matrix<double, Dynamic, Dynamic> , Lower> LPi_diag[Nhoriz]; // ***** [horizon]
    Matrix<double, _n, _n> LPi_diag_transp[_N];
    LLT<Matrix<double, Dynamic, Dynamic> , Lower> LRho_diag[Nhoriz]; // ***** [horizon]
    Matrix<double, _n, _n> LRho_diag_transp[_N];
    Matrix<Type, _m, _n> LSigma_offDiag[_N]; // last element used only if _pos_omega == _N
    Matrix<Type, _n, _m> LSigma_offDiag_transp[_N]; // last element used only if _pos_omega == _N
    Matrix<Type, _m, _n> LLambda0; // used for LLT decomposition if _pos_omega != _N
    Matrix<Type, _n, _m> LLambda0_transp;
    Matrix<Type, _m, _m> LLambda1; // used for LLT decomposition if _pos_omega != _N
    Matrix<Type, _m, _m> LLambda1_transp;

    // Variables needed to compute Y = C'*Phi_tilde * C;
    Matrix<Type, _m, _n> U[3 + (_N - 1) * 3]; // L*U = C'; last element only used in case _pos_omega == _N
    Matrix<Type, _n, _n> U_bar[2 + (_N - 1) * 5]; // L*U = C'
    Matrix<Type, _m, _n> UO[3]; // off-diag elements for _pos_omega != _N
    Matrix<Type, _m, _n> X[3 + (_N - 1) * 3]; // L'*X = U
    Matrix<Type, _n, _n> X_bar[2 + (_N - 1) * 5]; // L'*X = U
    Matrix<Type, _m, _n> XO[3]; // off-diag elements for _pos_omega != _N

    Matrix<Type, 2 * _N * _n, 1> beta; // beta = r_C + C*(Phi\r_d_bar), where r_d_bar = -r_H - P'*(lambda.*r_P./t - r_T./t);

    // pointers might be dangerous b/c no control over length
    const Matrix<Type, _n, 1> *x_star; // at each time step, we would like to track x_star[]
    Matrix<Type, _m, 1> u_star[_N]; // u_star[] is the least squares problem to track x_star[], but might not be feasible

    // ---------- private methods ----------
    void compRQ(); // computes the cost vectors r and q
    void compInitPoints(); // compute suitable initial starting points (z0, lambda0, nu0, slack0)
    void compResiduals(); // computes all 4 residua
    void compPhi();
    void compPhi_tilde(); // compute Cholesky decomposition
    //void compY(); // Schur complement
    inline void compY()
    {
      // computation of Y is done in three steps: Y = C * X
      // 1. L*U = C'
      // 2. L'*X = U --> L*L'*X = Phi*X = C', i.e. X = Phi_tilde * C'
      // 3. C*X = C*Phi_tilde*C' = Y

      Matrix<Type, _n, _n> eye;
      eye.setIdentity();

      // 1. Compute elements of Matrix U
      // treat the first 2 U_bar and first 3 U specially
      U[0] = LOmicron_diag[0].matrixLLT().triangularView<Lower> ().solve(-Bm_transp);
      U_bar[0] = LPi_diag[0].matrixLLT().triangularView<Lower> ().solve(eye);

      U[1] = LOmicron_diag[0].matrixLLT().triangularView<Lower> ().solve(-B_transp);
      U_bar[1] = LRho_diag[0].matrixLLT().triangularView<Lower> ().solve(eye);
      U[2] = LOmicron_diag[1].matrixLLT().triangularView<Lower> ().solve( /*zero*/-LSigma_offDiag[0] * U_bar[1]);

      // remaining U_bar and U in middle have structures
      for (int i = 1; i <= _N - 2; i++)
      {
        U_bar[2 + (i - 1) * 5] = LPi_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(-Am_tilde_transp);
        U_bar[3 + (i - 1) * 5] = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(-Bm_bar_transp);
        U[3 + (i - 1) * 3] = LOmicron_diag[i].matrixLLT().triangularView<Lower> ().solve(
                                                                                         -Bm_transp - LSigma_offDiag[i
                                                                                             - 1] * U_bar[3 + (i - 1)
                                                                                             * 5]);
        U_bar[4 + (i - 1) * 5] = LPi_diag[i].matrixLLT().triangularView<Lower> ().solve(eye);

        U_bar[5 + (i - 1) * 5] = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(-A_bar_transp);
        U[4 + (i - 1) * 3] = LOmicron_diag[i].matrixLLT().triangularView<Lower> ().solve(
                                                                                         -B_transp - LSigma_offDiag[i
                                                                                             - 1] * U_bar[5 + (i - 1)
                                                                                             * 5]);
        U_bar[6 + (i - 1) * 5] = LRho_diag[i].matrixLLT().triangularView<Lower> ().solve(eye);
        U[5 + (i - 1) * 3] = LOmicron_diag[i + 1].matrixLLT().triangularView<Lower> ().solve(
        /*zero*/-LSigma_offDiag[i] * U_bar[6 + (i - 1) * 5]);
      }

      U_bar[2 + (_N - 2) * 5] = LPi_diag[_N - 2].matrixLLT().triangularView<Lower> ().solve(-Am_tilde_transp);
      U_bar[3 + (_N - 2) * 5] = LRho_diag[_N - 2].matrixLLT().triangularView<Lower> ().solve(-Bm_bar_transp);
      U[3 + (_N - 2) * 3] = LOmicron_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                                             -Bm_transp
                                                                                                 - LSigma_offDiag[_N
                                                                                                     - 2] * U_bar[3
                                                                                                     + (_N - 2) * 5]);
      U_bar[4 + (_N - 2) * 5] = LPi_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(eye);

      U_bar[5 + (_N - 2) * 5] = LRho_diag[_N - 2].matrixLLT().triangularView<Lower> ().solve(-A_bar_transp);
      U[4 + (_N - 2) * 3] = LOmicron_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                                             -B_transp
                                                                                                 - LSigma_offDiag[_N
                                                                                                     - 2] * U_bar[5
                                                                                                     + (_N - 2) * 5]);
      U_bar[6 + (_N - 2) * 5] = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(eye);

      if (_N == _pos_omega)
        U[5 + (_N - 2) * 3] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
        /*zero*/-LSigma_offDiag[_N - 1] * U_bar[6 + (_N - 2) * 5]);
      else
      {
        UO[0] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
        /*zero*/-LLambda0 * U_bar[1 + (_pos_omega - 1) * 5] - LLambda1 * U[2 + (_pos_omega - 1) * 3]);
        UO[1] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
        /*zero*/-LLambda0 * U_bar[3 + (_pos_omega - 1) * 5] - LLambda1 * U[3 + (_pos_omega - 1) * 3]);
        UO[2] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
        /*zero*/-LLambda0 * U_bar[5 + (_pos_omega - 1) * 5] - LLambda1 * U[4 + (_pos_omega - 1) * 3]);
      }

      /*
       cout << "U[0]" << endl << U[0] << endl << endl;
       cout << "U_bar[0]" << endl << U_bar[0] << endl << endl;
       cout << "U[1]" << endl << U[1] << endl << endl;
       cout << "U_bar[1]" << endl << U_bar[1] << endl << endl;
       cout << "U[2]" << endl << U[2] << endl << endl;
       for (int i=1; i<= _N-2; i++)
       {
       cout << "U_bar[" << 2+(i-1)*5 << "]" << endl << U_bar[2+(i-1)*5] << endl << endl ;
       cout << "U_bar[" << 3+(i-1)*5 << "]" << endl << U_bar[3+(i-1)*5] << endl << endl ;
       cout << "U[" << 3+(i-1)*3 << "]" << endl << U[3+(i-1)*3] << endl << endl;
       cout << "U_bar[" << 4+(i-1)*5 << "]" << endl << U_bar[4+(i-1)*5] << endl << endl;

       cout << "U_bar[" << 5+(i-1)*5 << "]" << endl << U_bar[5+(i-1)*5] << endl << endl ;
       cout << "U[" << 4+(i-1)*3 << "]" << endl << U[4+(i-1)*3] << endl << endl ;
       cout << "U_bar[" << 6+(i-1)*5 << "]" << endl << U_bar[6+(i-1)*5] << endl << endl ;
       cout << "U[" << 5+(i-1)*3 << "]" << endl << U[5+(i-1)*3] << endl << endl ;
       }
       cout << "U_bar[" << 2+(_N-2)*5 << "]" << endl << U_bar[2+(_N-2)*5] << endl << endl ;
       cout << "U_bar[" << 3+(_N-2)*5 << "]" << endl << U_bar[3+(_N-2)*5] << endl << endl ;
       cout << "U[" << 3+(_N-2)*3 << "]" << endl << U[3+(_N-2)*3] << endl << endl;
       cout << "U_bar[" << 4+(_N-2)*5 << "]" << endl << U_bar[4+(_N-2)*5] << endl << endl;

       cout << "U_bar[" << 5+(_N-2)*5 << "]" << endl << U_bar[5+(_N-2)*5] << endl << endl ;
       cout << "U[" << 4+(_N-2)*3 << "]" << endl << U[4+(_N-2)*3] << endl << endl ;
       cout << "U_bar[" << 6+(_N-2)*5 << "]" << endl << U_bar[6+(_N-2)*5] << endl << endl ;

       if (_N == _pos_omega)
       {
       cout << "U[" << 5+(_N-2)*3 << "]" << endl << U[5+(_N-2)*3] << endl << endl ;
       }
       else
       {
       cout << "UO[0]" << endl << UO[0] << endl << endl;
       cout << "UO[1]" << endl << UO[1] << endl << endl;
       cout << "UO[2]" << endl << UO[2] << endl << endl;
       }
       */

      // 2. Compute elements in Matrix X
      // treat the first 2 X_bar and first 3 X specially

      if (_pos_omega != _N)
      {
        XO[0] = LOmicron_diag_transp[_N].template triangularView<Upper> ().solve(UO[0]);
        XO[1] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[1]);
        XO[2] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[2]);
      }

      X[0] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[0]);
      X_bar[0] = LPi_diag_transp[0].template triangularView<Upper>().solve(U_bar[0]);

      if (_pos_omega == 1) // build col1

            {
              X[1] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[1]);
              X[2] = LOmicron_diag_transp[1].template triangularView<Upper>().solve(U[2] - LLambda1_transp*XO[0]);
              X_bar[1] = LRho_diag_transp[0].template triangularView<Upper>().solve(U_bar[1] - LSigma_offDiag_transp[0]*X[2] - LLambda0_transp*XO[0]);
            }
            else
            {
              X[1] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[1]);
              X[2] = LOmicron_diag_transp[1].template triangularView<Upper>().solve(U[2]);
              X_bar[1] = LRho_diag_transp[0].template triangularView<Upper>().solve(U_bar[1] - LSigma_offDiag_transp[0]*X[2]);
            }

            // remaining X_bar and X have structures
            // if(_N != _pos_omega), then the off-diagonal element incluences three columns, i.e col1, col2, col3
            for (int i = 1; i <= _N-2; i++)
            {
              if (i == _pos_omega) // col2

              {
                X_bar[2+(_pos_omega-1)*5] = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[2+(_pos_omega-1)*5]);
                X[3+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(U[3+(_pos_omega-1)*3]-LLambda1_transp*XO[1]);
                X_bar[3+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[3+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[3+(_pos_omega-1)*3] - LLambda0_transp*XO[1]);
                X_bar[4+(_pos_omega-1)*5] = LPi_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar[4+(_pos_omega-1)*5]);
              }
              else // standard

              {
                X_bar[2+(i-1)*5] = LPi_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[2+(i-1)*5]);
                X[3+(i-1)*3] = LOmicron_diag_transp[i].template triangularView<Upper>().solve(U[3+(i-1)*3]);
                X_bar[3+(i-1)*5] = LRho_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[3+(i-1)*5] - LSigma_offDiag_transp[i-1]*X[3+(i-1)*3]);
                X_bar[4+(i-1)*5] = LPi_diag_transp[i].template triangularView<Upper>().solve(U_bar[4+(i-1)*5]);
              }

              if (i == _pos_omega) // col3

              {
                X[4+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve( U[4+(_pos_omega-1)*3] - LLambda1_transp*XO[2]);
                X_bar[5+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[5+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[4+(_pos_omega-1)*3] - LLambda0_transp*XO[2]);
                X[5+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega+1].template triangularView<Upper>().solve(U[5+(_pos_omega-1)*3]);
                X_bar[6+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar[6+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega]*X[5+(_pos_omega-1)*3]);
              }
              else if(i == _pos_omega-1) // col1

              {
                X[4+(_pos_omega-2)*3] = LOmicron_diag_transp[_pos_omega-1].template triangularView<Upper>().solve( U[4+(_pos_omega-2)*3] );
                X_bar[5+(_pos_omega-2)*5] = LRho_diag_transp[_pos_omega-2].template triangularView<Upper>().solve(U_bar[5+(_pos_omega-2)*5] - LSigma_offDiag_transp[_pos_omega-2]*X[4+(_pos_omega-2)*3]);
                X[5+(_pos_omega-2)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(U[5+(_pos_omega-2)*3] - LLambda1_transp*XO[0]);
                X_bar[6+(_pos_omega-2)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[6+(_pos_omega-2)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[5+(_pos_omega-2)*3] - LLambda0_transp*XO[0]);
              }
              else // if the off-diag element in P*z<h has no influence

              {
                X[4+(i-1)*3] = LOmicron_diag_transp[i].template triangularView<Upper>().solve( U[4+(i-1)*3] );
                X_bar[5+(i-1)*5] = LRho_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[5+(i-1)*5] - LSigma_offDiag_transp[i-1]*X[4+(i-1)*3]);
                X[5+(i-1)*3] = LOmicron_diag_transp[i+1].template triangularView<Upper>().solve(U[5+(i-1)*3]);
                X_bar[6+(i-1)*5] = LRho_diag_transp[i].template triangularView<Upper>().solve(U_bar[6+(i-1)*5] - LSigma_offDiag_transp[i]*X[5+(i-1)*3]);
              }
            }

            // compute last two columns
            if (_pos_omega == _N)
            {
              X_bar[2+(_N-2)*5] = LPi_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[2+(_N-2)*5]);
              X[3+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve(U[3+(_N-2)*3]);
              X_bar[3+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[3+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[3+(_N-2)*3]);
              X_bar[4+(_N-2)*5] = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[4+(_N-2)*5]);

              X[4+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve( U[4+(_N-2)*3] );
              X_bar[5+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[5+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[4+(_N-2)*3]);
              X[5+(_N-2)*3] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(U[5+(_N-2)*3]);
              X_bar[6+(_N-2)*5] = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[6+(_N-2)*5] - LSigma_offDiag_transp[_N-1]*X[5+(_N-2)*3]);
            }
            else if(_pos_omega == _N-1) // compute col2 and col3

            {
              X_bar[2+(_N-2)*5] = LPi_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[2+(_N-2)*5]);
              X[3+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve(U[3+(_N-2)*3] - LLambda1_transp*XO[1]);
              X_bar[3+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[3+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[3+(_N-2)*3] - LLambda0_transp*XO[1]);
              X_bar[4+(_N-2)*5] = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[4+(_N-2)*5]);

              X[4+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve( U[4+(_N-2)*3] - LLambda1_transp*XO[2]);
              X_bar[5+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[5+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[4+(_N-2)*3] - LLambda0_transp*XO[2]);
              X_bar[6+(_N-2)*5] = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[6+(_N-2)*5]);
            }
            else // standard, no influence by off-diag element in P

            {
              X_bar[2+(_N-2)*5] = LPi_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[2+(_N-2)*5]);
              X[3+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve(U[3+(_N-2)*3]);
              X_bar[3+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[3+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[3+(_N-2)*3]);
              X_bar[4+(_N-2)*5] = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[4+(_N-2)*5]);

              X[4+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve( U[4+(_N-2)*3]);
              X_bar[5+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[5+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[4+(_N-2)*3]);
              X_bar[6+(_N-2)*5] = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[6+(_N-2)*5]);
            }

            /*
             cout << "X[0]" << endl << X[0] << endl << endl;
             cout << "X_bar[0]" << endl << X_bar[0] << endl << endl;
             cout << "X[1]" << endl << X[1] << endl << endl;
             cout << "X_bar[1]" << endl << X_bar[1] << endl << endl;
             cout << "X[2]" << endl << X[2] << endl << endl;
             for (int i=1; i<= _N-2; i++)
             {
             cout << "X_bar[" << 2+(i-1)*5 << "]" << endl << X_bar[2+(i-1)*5] << endl << endl ;
             cout << "X_bar[" << 3+(i-1)*5 << "]" << endl << X_bar[3+(i-1)*5] << endl << endl ;
             cout << "X[" << 3+(i-1)*3 << "]" << endl << X[3+(i-1)*3] << endl << endl;
             cout << "X_bar[" << 4+(i-1)*5 << "]" << endl << X_bar[4+(i-1)*5] << endl << endl;

             cout << "X_bar[" << 5+(i-1)*5 << "]" << endl << X_bar[5+(i-1)*5] << endl << endl ;
             cout << "X[" << 4+(i-1)*3 << "]" << endl << X[4+(i-1)*3] << endl << endl ;
             cout << "X_bar[" << 6+(i-1)*5 << "]" << endl << X_bar[6+(i-1)*5] << endl << endl ;
             cout << "X[" << 5+(i-1)*3 << "]" << endl << X[5+(i-1)*3] << endl << endl ;
             }
             cout << "X_bar[2+(_N-2)*5]" << endl << X_bar[2+(_N-2)*5] << endl << endl ;
             cout << "X_bar[3+(_N-2)*5]" << endl << X_bar[3+(_N-2)*5] << endl << endl ;
             cout << "X[3+(_N-2)*3]" << endl << X[3+(_N-2)*3] << endl << endl;
             cout << "X_bar[4+(_N-2)*5]" << endl << X_bar[4+(_N-2)*5] << endl << endl;

             cout << "X_bar[5+(_N-2)*5]" << endl << X_bar[5+(_N-2)*5] << endl << endl ;
             cout << "X[4+(_N-2)*3]" << endl << X[4+(_N-2)*3] << endl << endl ;
             cout << "X_bar[6+(_N-2)*5]" << endl << X_bar[6+(_N-2)*5] << endl << endl ;
             if (_pos_omega == _N)
             {
             cout << "X[5+(_N-2)*3]" << endl << X[5+(_N-2)*3] << endl << endl ;
             }
             else
             {
             cout << "XO[0]" << endl << XO[0] << endl << endl;
             cout << "XO[1]" << endl << XO[1] << endl << endl;
             cout << "XO[2]" << endl << XO[2] << endl << endl;
             }
             */

            // 3. Compute Y = C*X
            // compute first three Y separately
            Y[0][0] = -Bm*X[0] + X_bar[0];

            Y[0][1] = -Bm*X[1];
            Y[1][1] = -B*X[1] + X_bar[1];

            // compute rest by filling Y column by column; done for 2 neighboring rows
            // we fill Y column by column, treating the first two columns specially
            for (int i=1; i <= _N-1; i++)
            {
              Y[0][2*(i-1)+2] = X_bar[2+(i-1)*5];
              Y[1][2*(i-1)+2] = X_bar[3+(i-1)*5];
              Y[2][2*(i-1)+2] = -Am_tilde*X_bar[2+(i-1)*5] - Bm_bar*X_bar[3+(i-1)*5] - Bm*X[3+(i-1)*3] + X_bar[4+(i-1)*5];

              Y[0][2*(i-1)+3] = X_bar[5+(i-1)*5];
              Y[1][2*(i-1)+3] = -Bm_bar*X_bar[5+(i-1)*5] - Bm*X[4+(i-1)*3];
              Y[2][2*(i-1)+3] = -A_bar*X_bar[5+(i-1)*5] - B*X[4+(i-1)*3] + X_bar[6+(i-1)*5];
            }

            /*
             cout << "Y[0][0]" << endl << Y[0][0] << endl << endl;
             cout << "Y[0][1]" << endl << Y[0][1] << endl << endl;
             cout << "Y[1][1]" << endl << Y[1][1] << endl << endl;
             for (int i=1; i <= _N-1; i++)
             {
             cout << "Y[0][2*(i-1)+2]" << endl << Y[0][2*(i-1)+2] << endl << endl;
             cout << "Y[1][2*(i-1)+2]" << endl << Y[1][2*(i-1)+2] << endl << endl;
             cout << "Y[2][2*(i-1)+2]" << endl << Y[2][2*(i-1)+2] << endl << endl;

             cout << "Y[0][2*(i-1)+3]" << endl << Y[0][2*(i-1)+3] << endl << endl;
             cout << "Y[1][2*(i-1)+3]" << endl << Y[1][2*(i-1)+3] << endl << endl;
             cout << "Y[2][2*(i-1)+3]" << endl << Y[2][2*(i-1)+3] << endl << endl;
             }
             */
          }

          void compBeta(); // beta = -r_p + C*Phi_tilde*r_d;
          void compL(); // L*L' = Y
          void compDnu();
          void compDz();
          void compDlambda(); // dlambda = (P*dz - r_P + r_T./lambda) .* lambda ./ t
          void compAlpha();
          double getTimer(); // method returning the computational time

        public:
          LBmpcTP(); // default constructor
          // standard constructor
          LBmpcTP(const char fileName[], bool verbose_arg);

          // "step" computes and returns status code
          int step(const Matrix<Type, _n, _n> &Lm_arg, const Matrix<Type, _n, _m> &Mm_arg, const Matrix<Type, _n, 1> &tm_arg,
              const Matrix<Type, _n, 1> &x_hat_arg, const Matrix<Type, _n, 1> x_star_arg[]);

          Matrix<Type, _m, 1> u_opt;
          int n_iter_last; // number of iterations in previous step
          //~LBmpcTP();	// destructor

          void get_A(Matrix<Type, _n, _n>& Aout)
          {
            Aout = A;
          }

          void get_B(Matrix<Type, _n, _m>& Bout)
          {
            Bout = B;
          }

          void get_s(Matrix<Type, _n, 1>& sout)
          {
            sout = s;
          }
        };

        //  ==================== Implementation of Methods ==================
        // default constructor
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::LBmpcTP()
  {
    cerr << "LBmpcTP() constructor was called with no arguments -- Please provide fileName" << endl;
  }

// constructor
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::LBmpcTP //
  (const char fileName[], bool verbose_arg) :
    n_iter_last(0)
  {
    // first starting z can be modified
    z.setConstant(0);
    nu.setConstant(1); // any nu good to initialize
    lambda.setConstant(1);
    slack.setConstant(1);
    verbose = verbose_arg;
    // srand ( time(NULL) );
    // nu.setRandom(); nu = 100*nu;
    // lambda.setRandom(); lambda = 100*lambda;
    // slack.setRandom(); slack = 100*slack;

    // -------------- read from binary file -----------------
    ifstream fin; // Definition input file object
    fin.open(fileName, ios::binary); //	open file
    if (!fin.is_open())
    {
      cerr << "File open error \n";
      return;
    }

    // read
    fin.read((char *)&n_iter, sizeof(int));
    fin.read((char *)&reg, sizeof(double));
    fin.read((char *)&resNorm_H, sizeof(double));
    fin.read((char *)&resNorm_C, sizeof(double));
    fin.read((char *)&resNorm_P, sizeof(double));
    fin.read((char *)&muNorm, sizeof(double));

    // read A
    for (int i = 0; i <= _n - 1; i++)
      for (int j = 0; j <= _n - 1; j++)
        fin.read((char *)&A(j, i), sizeof(double));

    // read B
    for (int i = 0; i <= _m - 1; i++) // #columns
      for (int j = 0; j <= _n - 1; j++) // #rows
        fin.read((char *)&B(j, i), sizeof(double));

    // read s
    for (int i = 0; i <= _n - 1; i++) // #columns
      fin.read((char *)&s(i, 0), sizeof(double));

    // read Q_tilde
    for (int i = 0; i <= _n - 1; i++) // #columns
      for (int j = 0; j <= _n - 1; j++) // #rows
        fin.read((char *)&Q_tilde(j, i), sizeof(double));

    // read Q_tilde_f
    for (int i = 0; i <= _n - 1; i++) // #columns
      for (int j = 0; j <= _n - 1; j++) // #rows
        fin.read((char *)&Q_tilde_f(j, i), sizeof(double));

    // read R
    for (int i = 0; i <= _m - 1; i++) // #columns
      for (int j = 0; j <= _m - 1; j++) // #rows
        fin.read((char *)&R(j, i), sizeof(double));

    // read Fx[_N]
    for (int k = 0; k <= _N - 1; k++)
      for (int i = 0; i <= _n - 1; i++) // #columns
        for (int j = 0; j <= _nSt - 1; j++) // #rows
          fin.read((char *)&Fx[k](j, i), sizeof(double));

    // read fx[_N]
    for (int k = 0; k <= _N - 1; k++)
      for (int i = 0; i <= _nSt - 1; i++) // #columns
        fin.read((char *)&fx[k](i, 0), sizeof(double));

    // read Fu[_N]
    for (int k = 0; k <= _N - 1; k++)
      for (int i = 0; i <= _m - 1; i++) // #columns
        for (int j = 0; j <= _nInp - 1; j++) // #rows
          fin.read((char *)&Fu[k](j, i), sizeof(double));

    // read fu[_N]
    for (int k = 0; k <= _N - 1; k++)
      for (int i = 0; i <= _nInp - 1; i++) // #columns
        fin.read((char *)&fu[k](i, 0), sizeof(double));

    // read F_xTheta
    for (int i = 0; i <= _n - 1; i++) // #columns
      for (int j = 0; j <= _nF_xTheta - 1; j++) // #rows
        fin.read((char *)&F_xTheta(j, i), sizeof(double));

    // read F_theta
    for (int i = 0; i <= _m - 1; i++) // #columns
      for (int j = 0; j <= _nF_xTheta - 1; j++) // #rows
        fin.read((char *)&F_theta(j, i), sizeof(double));

    // read f_xTheta
    for (int i = 0; i <= _nF_xTheta - 1; i++) // #columns
      fin.read((char *)&f_xTheta(i, 0), sizeof(double));

    // read K_arg
    for (int i = 0; i <= _n - 1; i++) // #columns
      for (int j = 0; j <= _m - 1; j++) // #rows
        fin.read((char *)&K(j, i), sizeof(double));

    fin.close(); // close file

    damp = 0.95;

    offset = _n + _n + _m;
    A_transp = A.transpose();
    B_transp = B.transpose();
    K_transp = K.transpose();

    for (int i = 0; i < _N; i++)
    {
      Fx_transp[i] = Fx[i].transpose();
      Fu_transp[i] = Fu[i].transpose();
      Fu_bar[i] = Fu[i] * K; // the first one is NOT used
      Fu_bar_transp[i] = Fu_bar[i].transpose();
    }

    F_xTheta_transp = F_xTheta.transpose();
    F_theta_transp = F_theta.transpose();

    // do some preliminary Matrix calculations
    S = K.transpose() * R;
    S_transp = S.transpose();
    Q_bar = S * K;
    A_bar = A + B * K;
    A_bar_transp = A_bar.transpose();

    /*
     // test output
     cout << "n_iter: " << n_iter << endl << endl;
     cout << "reg: " << reg << endl << endl;
     cout << "resNorm_H: " << resNorm_H << endl << endl;
     cout << "resNorm_C: " << resNorm_C << endl << endl;
     cout << "resNorm_P: " << resNorm_P << endl << endl;
     cout << "muNorm: " << muNorm << endl << endl;

     cout << "A: " << endl << A << endl << endl;
     cout << "B: " << endl << B << endl << endl;
     cout << "s:" << endl << s << endl << endl;
     cout << "Q_tilde: " << endl << Q_tilde << endl << endl;
     cout << "Q_tilde_f: " << endl << Q_tilde_f << endl << endl;
     cout << "R: " << endl << R << endl << endl;

     for (int i = 0; i<= _N-1; i++)
     {
     cout << "Fx[" << i << "]: " << endl << Fx[i] << endl << endl;
     }
     for (int i = 0; i <= _N-1; i++)
     {
     cout << "fx[" << i << "]: " << endl << fx[i] << endl << endl;
     }
     for (int i = 0; i <= _N-1; i++)
     {
     cout << "Fu[" << i << "]: " << endl << Fu[0] << endl << endl;
     }
     for (int i = 0; i <= _N-1; i++)
     {
     cout << "fu[" << i << "]: " << endl << fu[i] << endl << endl;
     }
     cout << "F_xTheta: " << endl << F_xTheta << endl << endl;
     cout << "F_theta: " << endl << F_theta << endl << endl;
     cout << "f_xTheta: " << endl << f_xTheta << endl << endl;
     cout << "K:" << endl << K << endl << endl;
     */
  }

// step function returns optimal input
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  int LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::step(const Matrix<Type, _n, _n> &Lm_arg,
                                                                           const Matrix<Type, _n, _m> &Mm_arg,
                                                                           const Matrix<Type, _n, 1> &tm_arg,
                                                                           const Matrix<Type, _n, 1> &x_hat_arg,
                                                                           const Matrix<Type, _n, 1> x_star_arg[])
  {
    if (verbose)
      cout << "*********** START OPTIMIZATION ***********" << endl;

    // initialization
    x_hat = &x_hat_arg;

    x_star = x_star_arg;
    tm_tilde = tm_arg + s;
    Am_tilde = A + Lm_arg;
    Am_tilde_transp = Am_tilde.transpose();
    Bm = B + Mm_arg;
    Bm_transp = Bm.transpose();
    Bm_bar = Bm * K;
    Bm_bar_transp = Bm_bar.transpose();

    //update z from previous z, which was assumed to be optimal or created during instantiation of object
    z.template segment<(_N - 1) * (_m + _n + _n)> (0) = z.template segment<(_N - 1) * (_m + _n + _n)> (offset);
    // ONLY TEMPORARILY
    // z << -2 ,    3  ,   3  ,   1 ,   -2 ,    2 ,   -3 ,    0  ,  -1  ,   0  ,   5 ,    3  ,   5  ,  -3  ,   0  ,  -4  ,   3  ,   1  ,   4  ,   5  ,   4 ,   -1,
    // -5  ,   0  ,  -3 ,   -3  ,  -2  ,  -4  ,   2  ,   2   ,  0  ,  -2  ,   3  ,   1 ,    5  ,   4  ,  -1  ,   0   , -2   ,  1  ,   3  ,   2 ,   -4   ,  3,
    // -5   , -1   ,  2   ,  3  ,  -1   ,  2;
    // srand ( time(NULL) );
    // z.setRandom();
    // z = 100*z;
    compRQ(); // compute u_star, x_star -> cost matrices
    compInitPoints(); // computes more suitable initial points, also for infeasible start
    // but makes warm start less useful

    int itNewton = 0;
    double sigma;
    compResiduals(); // compute r_H, r_C, r_P, r_T
    double mu_pd = lambda.dot(slack) / (_N * (_nInp + _nSt) + _nF_xTheta); // mu of primal-dual IPM
    double mu_pd_p; // mu_p of for predictor step in primal-dual IPM
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> ones;
    bool condition;
    do
    {
      itNewton++;
      if (itNewton >= n_iter)
      {
        n_iter_last = itNewton - 1;
        if (verbose)
          cerr << "-------- more than " << n_iter << " Newton steps required" << endl;
        if (r_H.norm() > resNorm_H)
        {
          if (verbose)
            cout << "norm(r_H) > " << resNorm_H << endl;
          return 1;
        }
        else if (r_C.norm() > resNorm_C)
        {
          if (verbose)
            cout << "norm(r_C) > " << resNorm_C << endl;
          return 2;
        }
        else if (r_P.norm() > resNorm_P)
        {
          if (verbose)
            cout << "norm(r_P) > " << resNorm_P << endl;
          return 3;
        }
        else if (mu_pd > muNorm)
        {
          if (verbose)
            cout << "mu > " << muNorm << endl;
          return 4;
        }
        else
          return 6;
        break;
      }
      // compute dz, dnu, dlambda, dslack
      compPhi(); // Phi = 2*H + kappa*P'*diag(1./slack)*diag(lambda)*P
      compPhi_tilde(); // Phi_tilde = Phi^{-1}; computes Chol-Dec
      compY(); // Y = C*Phi_tilde*C'
      compBeta(); // beta = r_C + C*(Phi\r_d_bar), where r_d_bar = -r_H - P'*(lambda.*r_P./t - r_T./t);
      compL(); // L*L' = Y
      compDnu(); // Y * dnu = -beta
      compDz(); // Phi*dz = -r_d_bar - C'*dnu
      compDlambda(); // dlambda = (P*dz - r_P + r_T./lambda) .* lambda ./ t
      dslack = (r_T - slack.cwiseProduct(dlambda)).cwiseQuotient(lambda); // compute dslack = (r_T - slack.*dlambda) ./ lambda

      compAlpha(); // computes steplength for predictor
      mu_pd_p = (lambda + alpha * dlambda).dot(slack + alpha * dslack) / (_N * (_nSt + _nInp) + _nF_xTheta);
      sigma = pow(mu_pd_p / mu_pd, 3);
      r_T = r_T - dlambda.cwiseProduct(dslack) + ones.setConstant(sigma * mu_pd);
      compBeta();
      compDnu(); // Y * dnu = -beta
      compDz(); // Phi*dz = -r_d_bar - C'*dnu

      compDlambda(); // dlambda = (P*dz - r_P + r_T./lambda) .* lambda ./ t
      dslack = (r_T - slack.cwiseProduct(dlambda)).cwiseQuotient(lambda); // compute dslack = (r_T - slack.*dlambda) ./ lambda
      compAlpha();

      // update z, nu, lambda, slack
      for (int i = 0; i <= _N * (_m + _n + _n) + _m - 1; i++)
      {
        if (isnan(dz[i]))
        {
          if (verbose)
            cout << "NAN detected during iteration number " << itNewton << endl;
          u_opt = K * (*x_hat) + z.template segment<_m> (0);
          return 5;
        }
      }
      z = z + damp * alpha * dz;

      nu = nu + damp * alpha * dnu;
      lambda = lambda + damp * alpha * dlambda;
      slack = slack + damp * alpha * dslack;

      compResiduals(); // compute r_H, r_C, r_P, r_T
      mu_pd = lambda.dot(slack) / (_N * (_nInp + _nSt) + _nF_xTheta);
      condition = (r_H.norm() > resNorm_H) || (r_C.norm() > resNorm_C) || (r_P.norm() > resNorm_P) || (mu_pd > muNorm);
    } while (condition);

    if (verbose)
    {
      cout << " =====> computed optimal z_vector:" << endl << setprecision(30) << z << endl << endl;
      cout << "number of Newton iterations required: " << itNewton << endl << endl;
    }
    u_opt = K * (*x_hat) + z.template segment<_m> (0);
    n_iter_last = itNewton;
    return 0;
  }

// ------------ function computes residua r_H, r_C, r_P, r_T -------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compResiduals()
  {
    // 1. compute r_H = -2*H*z - g - P'*lambda - C'*nu;
    // handle first case separately
    r_H.template segment<_m> (0) = 2 * R * z.template segment<_m> (0) + (r_vec[0] + 2 * S_transp * (*x_hat))
        + Fu_transp[0] * lambda.template segment<_nInp> (0) - Bm_transp * nu.template segment<_n> (0) - B_transp
        * nu.template segment<_n> (_n);

    // handle the cases in the middle, without the end, three block to deal with in each round
    int offset1 = 2 * _n; // offset required in nu for C'*nu
    for (int i = 1; i <= _N - 1; i++)
    {
      r_H.template segment<_n> (_m + (i - 1) * offset) = 2 * Q_tilde * z.template segment<_n> (_m + (i - 1) * offset)
          + q_tilde_vec[i - 1] + nu.template segment<_n> ((i - 1) * offset1) - Am_tilde_transp
          * nu.template segment<_n> ((i - 1) * offset1 + _n + _n);

      if (i != _pos_omega)
      {
        r_H.template segment<_n> (_m + (i - 1) * offset + _n) = 2 * Q_bar * z.template segment<_n> (
                                                                                                    _m + (i - 1)
                                                                                                        * offset + _n)
            + 2 * S * z.template segment<_m> (_m + (i - 1) * offset + _n + _n) + q_bar_vec[i] + Fx_transp[i - 1]
            * lambda.template segment<_nSt> (_nInp + (i - 1) * (_nInp + _nSt)) + Fu_bar_transp[i]
            * lambda.template segment<_nInp> (i * (_nInp + _nSt)) + nu.template segment<_n> ((i - 1) * offset1 + _n)
            - Bm_bar_transp * nu.template segment<_n> ((i - 1) * offset1 + _n + _n) - A_bar_transp
            * nu.template segment<_n> ((i - 1) * offset1 + _n + _n + _n);
      }
      else // must add the additional term: F_xTheta'*d(.)
      {
        r_H.template segment<_n> (_m + (_pos_omega - 1) * offset + _n) = 2 * Q_bar
            * z.template segment<_n> (_m + (_pos_omega - 1) * offset + _n) + 2 * S
            * z.template segment<_m> (_m + (_pos_omega - 1) * offset + _n + _n) + q_bar_vec[i] + Fx_transp[_pos_omega
            - 1] * lambda.template segment<_nSt> (_nInp + (i - 1) * (_nInp + _nSt)) + Fu_bar_transp[_pos_omega]
            * lambda.template segment<_nInp> (i * (_nInp + _nSt)) + nu.template segment<_n> (
                                                                                             (_pos_omega - 1) * offset1
                                                                                                 + _n) - Bm_bar_transp
            * nu.template segment<_n> ((_pos_omega - 1) * offset1 + _n + _n) - A_bar_transp
            * nu.template segment<_n> ((_pos_omega - 1) * offset1 + _n + _n + _n) + F_xTheta_transp
            * lambda.template segment<_nF_xTheta> (_N * (_nSt + _nInp)); // used to be: num_constr - _nF_xTheta
      }

      r_H.template segment<_m> (_m + (i - 1) * offset + _n + _n) = 2 * S_transp * z.template segment<_n> (
                                                                                                          _m + (i - 1)
                                                                                                              * offset
                                                                                                              + _n) + 2
          * R * z.template segment<_m> (_m + (i - 1) * offset + _n + _n) + r_vec[i] + Fu_transp[i]
          * lambda.template segment<_nInp> (i * (_nInp + _nSt)) - Bm_transp * nu.template segment<_n> (
                                                                                                       (i - 1)
                                                                                                           * offset1
                                                                                                           + _n + _n)
          - B_transp * nu.template segment<_n> ((i - 1) * offset1 + _n + _n + _n);
    }
    r_H.template segment<_n> (_m + (_N - 1) * offset) = 2 * Q_tilde_f * z.template segment<_n> (_m + (_N - 1) * offset)
        + q_tilde_vec[_N - 1] + nu.template segment<_n> ((_N - 1) * offset1);

    if (_pos_omega == _N)
    {
      r_H.template segment<_n> (_m + (_N - 1) * offset + _n) = Fx_transp[_N - 1]
          * lambda.template segment<_nSt> (_nInp + (_N - 1) * (_nInp + _nSt)) + F_xTheta_transp
          * lambda.template segment<_nF_xTheta> (_N * (_nInp + _nSt)) + nu.template segment<_n> (
                                                                                                 (_N - 1) * offset1
                                                                                                     + _n);
    }
    else //standard
    {
      r_H.template segment<_n> (_m + (_N - 1) * offset + _n) = Fx_transp[_N - 1]
          * lambda.template segment<_nSt> (_nInp + (_N - 1) * (_nInp + _nSt)) + nu.template segment<_n> (
                                                                                                         (_N - 1)
                                                                                                             * offset1
                                                                                                             + _n);
    }
    r_H.template segment<_m> (_m + (_N - 1) * offset + _n + _n) = F_theta_transp
        * lambda.template segment<_nF_xTheta> (_N * (_nSt + _nInp));
    r_H = -r_H;
    // cout << setprecision(30) << "r_H" << endl << r_H << endl << endl;

    // 2. r_C = b - C*z;
    r_C.template segment<_n> (0) = -Bm * z.template segment<_m> (0) + z.template segment<_n> (_m) - (Am_tilde + Bm_bar)
        * (*x_hat) - tm_tilde;
    r_C.template segment<_n> (_n) = -B * z.template segment<_m> (0) + z.template segment<_n> (_m + _n) - (A_bar
        * (*x_hat) + s);

    // deal with the rest, class variable: offset = _n + _n + _m;
    for (int i = 1; i <= _N - 1; i++)
    {
      r_C.template segment<_n> (2 * _n * i) = -Am_tilde * z.template segment<_n> ((i - 1) * offset + _m) - Bm_bar
          * z.template segment<_n> ((i - 1) * offset + _m + _n) - Bm * z.template segment<_m> (
                                                                                               (i - 1) * offset + _m
                                                                                                   + _n + _n)
          + z.template segment<_n> ((i - 1) * offset + _m + _n + _n + _m) - tm_tilde;
      r_C.template segment<_n> (2 * _n * i + _n) = -A_bar * z.template segment<_n> ((i - 1) * offset + _m + _n) - B
          * z.template segment<_m> ((i - 1) * offset + _m + _n + _n) + z.template segment<_n> (
                                                                                               (i - 1) * offset + _m
                                                                                                   + _n + _n + _m + _n)
          - s;
    }
    r_C = -r_C;
    // cout << setprecision(30) << "r_C" << endl << r_C << endl << endl;


    // 3. r_P = -slack + h - P*z;
    c_tmp = z.template segment<_m> (0);
    r_P.template segment<_nInp> (0) = -slack.template segment<_nInp> (0) + (fu[0] - Fu[0] * K * (*x_hat)) - (Fu[0]
        * c_tmp); // should be >0

    // general treatment in the middle, class variable: offset = _n + _n + _m
    for (int i = 1; i <= _N - 1; i++) // blocks in the middle
    { // compute (h - P*z)_i
      x_bar_tmp = z.template segment<_n> ((i - 1) * offset + _m + _n);
      c_tmp = z.template segment<_m> ((i - 1) * offset + _m + _n + _n);
      r_P.template segment<_nSt> (_nInp + (i - 1) * (_nInp + _nSt)) = -slack.template segment<_nSt> (
                                                                                                     _nInp + (i - 1)
                                                                                                         * (_nInp
                                                                                                             + _nSt))
          + fx[i - 1] - Fx[i - 1] * x_bar_tmp;
      r_P.template segment<_nInp> (i * (_nInp + _nSt)) = -slack.template segment<_nInp> (i * (_nInp + _nSt)) + fu[i]
          - (Fu[i] * K * x_bar_tmp + Fu[i] * c_tmp);
    }
    // special case for last blocks
    x_bar_tmp = z.template segment<_n> ((_N - 1) * offset + _m + _n);
    c_tmp = z.template segment<_m> ((_N - 1) * offset + _m + _n + _n);
    r_P.template segment<_nSt> (_nInp + (_N - 1) * (_nSt + _nInp))
        = -slack.template segment<_nSt> (_nInp + (_N - 1) * (_nSt + _nInp)) + fx[_N - 1] - (Fx[_N - 1] * x_bar_tmp);
    x_bar_tmp = z.template segment<_n> ((_pos_omega - 1) * offset + _m + _n);
    r_P.template segment<_nF_xTheta> (_N * (_nSt + _nInp)) = -slack.template segment<_nF_xTheta> (_N * (_nSt + _nInp))
        + f_xTheta - (F_xTheta * x_bar_tmp + F_theta * c_tmp);
    // cout << setprecision(30) << "r_P:" << endl << r_P << endl << endl;

    // 4. r_T = -slack .* lambda
    r_T = -slack.cwiseProduct(lambda);
    // cout << "r_T:" << endl << r_T << endl << endl;
  }

// ------------ function computes suitable lambda0 and slack0  --------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compInitPoints()
  {

    // 1. compute dlambda and dslack
    // 2. lambda0 = max{1, |lambda+dlambda|}
    // 3. slack0 = max{1, |slack+dslack|}
    // 4. z0 = z; nu0 = nu;

    Array<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> arr_ones;
    arr_ones.setConstant(1);

    compResiduals(); // compute r_H, r_C, r_P, r_T
    compPhi(); // Phi = 2*H + kappa*P'*diag(1./slack)*diag(lambda)*P
    compPhi_tilde(); // Phi_tilde = Phi^{-1}; computes Chol-Dec
    compY(); // Y = C*Phi_tilde*C'
    compBeta(); // beta = r_C + C*(Phi\r_d_bar), where r_d_bar = -r_H - P'*(lambda.*r_P./t - r_T./t);
    compL(); // L*L' = Y
    compDnu(); // Y * dnu = -beta
    compDz(); // Phi*dz = -r_d_bar - C'*dnu
    compDlambda(); // dlambda = (P*dz - r_P + r_T./lambda) .* lambda ./ t
    dslack = (r_T - slack.cwiseProduct(dlambda)).cwiseQuotient(lambda); // compute dslack = (r_T - slack.*dlambda) ./ lambda

    // new slack and lambda
    slack = arr_ones.max((slack + dslack).array().abs());
    lambda = arr_ones.max((lambda + dlambda).array().abs());

    // cout << setprecision(30) << "new slack: " << endl << slack << endl << endl;
    // cout << "new lambda: " << endl << lambda << endl << endl;
    // cout << "new z: " << endl << z << endl << endl;
    // cout << "new nu: " << endl << nu << endl << endl;
  }

// ------- computes Phi = 2*H + P' * diag(1./slack) * diag(lambda) * P ------------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compPhi()
  {
    DiagonalMatrix<Type, _nInp> d_diagU[_N];
    DiagonalMatrix<Type, _nSt> d_diagX[_N];
    DiagonalMatrix<Type, _nF_xTheta> d_diagTheta;
    for (int i = 0; i <= _N - 1; i++)
    {
      d_diagU[i].diagonal()
          = (lambda.template segment<_nInp> (i * (_nInp + _nSt))).cwiseQuotient(
                                                                                slack.template segment<_nInp> (
                                                                                                               i
                                                                                                                   * (_nInp
                                                                                                                       + _nSt)));
      d_diagX[i].diagonal()
          = (lambda.template segment<_nSt> (_nInp + i * (_nInp + _nSt))).cwiseQuotient(
                                                                                       slack.template segment<_nSt> (
                                                                                                                     _nInp
                                                                                                                         + i
                                                                                                                             * (_nInp
                                                                                                                                 + _nSt)));
    }
    d_diagTheta.diagonal()
        = (lambda.template segment<_nF_xTheta> (_N * (_nInp + _nSt))).cwiseQuotient(
                                                                                    slack.template segment<_nF_xTheta> (
                                                                                                                        _N
                                                                                                                            * (_nInp
                                                                                                                                + _nSt)));

    // --------------- compute elements of Phi, i.e. elements of Omicron, Pi = 2*Q_tilde, Rho, Sigma
    Matrix<Type, _m, _m> eyeM;
    eyeM.setIdentity();
    Matrix<Type, _n, _n> eyeN;
    eyeN.setIdentity();

    // special treatment at the beginning
    Omicron[0] = 2 * R + Fu_transp[0] * d_diagU[0] * Fu[0] + reg * eyeM;

    // do the rest by computing three block and three block
    for (int i = 1; i <= _N - 1; i++)
    {
      if (i != _pos_omega)
      {
        Rho[i - 1] = 2 * Q_bar + Fx_transp[i - 1] * d_diagX[i - 1] * Fx[i - 1] + Fu_bar_transp[i] * d_diagU[i]
            * Fu_bar[i] + reg * eyeN;
      }
      else // i == _pos_omega
      {
        Rho[_pos_omega - 1] = 2 * Q_bar + (Fx_transp[_pos_omega - 1] * d_diagX[_pos_omega - 1] * Fx[i - 1]
            + Fu_bar_transp[_pos_omega] * d_diagU[_pos_omega] * Fu_bar[_pos_omega] + F_xTheta_transp * d_diagTheta
            * F_xTheta) + reg * eyeN;
      }
      Sigma[i - 1] = 2 * S + (Fu_bar_transp[i] * d_diagU[i] * Fu[i]);
      Omicron[i] = 2 * R + (Fu_transp[i] * d_diagU[i] * Fu[i]) + reg * eyeM;
    }

    // special treatment for last block
    if (_pos_omega == _N)
    {
      Rho[_N - 1] = (Fx_transp[_N - 1] * d_diagX[_N - 1] * Fx[_N - 1] + F_xTheta_transp * d_diagTheta * F_xTheta) + reg
          * eyeN;
    }
    else // considered in loop above
    {
      Rho[_N - 1] = (Fx_transp[_N - 1] * d_diagX[_N - 1] * Fx[_N - 1]) + reg * eyeN;
    }
    Sigma[_N - 1] = (F_xTheta_transp * d_diagTheta * F_theta); // independent of _pos_omega, represents the off-diag matrix
    Omicron[_N] = (F_theta_transp * d_diagTheta * F_theta) + reg * eyeM;

    /*
     cout << "reg: " << reg << endl << endl;
     for (int i = 0; i <= _N-1; i++)
     {
     cout << "Omicron[" << i << "]" << endl << Omicron[i] << endl << endl;
     cout << "Rho[" << i << "]" << endl << Rho[i] << endl << endl;
     cout << "Sigma[" << i << "]" << endl << Sigma[i] << endl << endl;
     }
     cout << "Omicron[" << _N << "]" << endl << Omicron[_N] << endl << endl;
     */
  }

// ------- computes Phi_tilde = Phi^{-1} -----------------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compPhi_tilde()
  {
    Matrix<Type, _n, _n> eyeN;
    eyeN.setIdentity();

    // decompose Phi = L_Phi*L_Phi'
    LOmicron_diag[0].compute(Omicron[0]);
    LOmicron_diag_transp[0] = LOmicron_diag[0].matrixLLT().transpose();
    for (int i = 1; i <= _N - 1; i++)
    {
      LPi_diag[i - 1].compute(2 * Q_tilde);
      LPi_diag_transp[i - 1] = LPi_diag[i - 1].matrixLLT().transpose();
      LRho_diag[i - 1].compute(Rho[i - 1]);
      LRho_diag_transp[i - 1] = LRho_diag[i - 1].matrixLLT().transpose();
      LSigma_offDiag_transp[i - 1] = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(Sigma[i - 1]);
      LSigma_offDiag[i - 1] = LSigma_offDiag_transp[i - 1].transpose();
      LOmicron_diag[i].compute(Omicron[i] - LSigma_offDiag[i - 1] * LSigma_offDiag_transp[i - 1]);
      LOmicron_diag_transp[i] = LOmicron_diag[i].matrixLLT().transpose();

      if (i == _pos_omega)
      {
        LLambda0_transp = LRho_diag[_pos_omega - 1].matrixLLT().triangularView<Lower> ().solve(Sigma[_N - 1]);
        LLambda0 = LLambda0_transp.transpose();
        LLambda1_transp
            = LOmicron_diag[_pos_omega].matrixLLT().triangularView<Lower> ().solve(
                                                                                   Matrix<Type, _m, _m>::Zero()
                                                                                       - LSigma_offDiag[_pos_omega - 1]
                                                                                           * LLambda0_transp);
        LLambda1 = LLambda1_transp.transpose();
      }
    }

    LPi_diag[_N - 1].compute(2 * Q_tilde_f);
    LPi_diag_transp[_N - 1] = LPi_diag[_N - 1].matrixLLT().transpose();
    LRho_diag[_N - 1].compute(Rho[_N - 1]);
    LRho_diag_transp[_N - 1] = LRho_diag[_N - 1].matrixLLT().transpose();

    if (_N == _pos_omega)
    {
      LSigma_offDiag_transp[_N - 1] = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(Sigma[_N - 1]);
      LSigma_offDiag[_N - 1] = LSigma_offDiag_transp[_N - 1].transpose();
      LOmicron_diag[_N].compute(Omicron[_N] - LSigma_offDiag[_N - 1] * LSigma_offDiag_transp[_N - 1]);
      LOmicron_diag_transp[_N] = LOmicron_diag[_N].matrixLLT().transpose();
    }
    else // LSigma_offDiag[_N-1] is not used
    {
      LOmicron_diag[_N].compute(Omicron[_N] - LLambda0 * LLambda0_transp - LLambda1 * LLambda1_transp);
      LOmicron_diag_transp[_N] = LOmicron_diag[_N].matrixLLT().transpose();
    }

    /*
     cout << "reg: " << reg << endl << endl;
     for (int i = 0; i <= _N-2; i++)
     {
     cout << "LOmicron_diag[" << i << "]" << endl << LOmicron_diag[i].matrixLLT() << endl << endl;
     cout << "LPi_diag[" << i << "]" << endl << LPi_diag[i].matrixLLT() << endl << endl;
     cout << "LRho_diag[" << i << "]" << endl << LRho_diag[i].matrixLLT() << endl << endl;
     cout << "LSigma_offDiag[" << i << "]" << endl << LSigma_offDiag[i] << endl << endl;
     }
     cout << "LOmicron_diag[" << _N-1 << "]" << endl << LOmicron_diag[_N-1].matrixLLT() << endl << endl;
     cout << "LPi_diag[" << _N-1 << "]" << endl << LPi_diag[_N-1].matrixLLT() << endl << endl;
     cout << "LRho_diag[" << _N-1 << "]" << endl << LRho_diag[_N-1].matrixLLT() << endl << endl;
     if (_N == _pos_omega)
     {
     cout << "LSigma_offDiag[" << _N-1 << "]" << endl << LSigma_offDiag[_N-1] << endl << endl;
     }
     cout << "LOmicron_diag[" << _N << "]" << endl << LOmicron_diag[_N].matrixLLT() << endl << endl;
     if (_pos_omega != _N)
     {
     cout << "LLambda0" << endl << LLambda0 << endl << endl;
     cout << "LLambda1" << endl << LLambda1 << endl << endl;
     }
     */

  }

// ------- computes Y = C*Phi_tilde*C' ------------------------
// ------------ works -----------------------------------------
//template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
//  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compY()
//  {
//    // computation of Y is done in three steps: Y = C * X
//    // 1. L*U = C'
//    // 2. L'*X = U --> L*L'*X = Phi*X = C', i.e. X = Phi_tilde * C'
//    // 3. C*X = C*Phi_tilde*C' = Y
//
//    Matrix<Type, _n, _n> eye;
//    eye.setIdentity();
//
//    // 1. Compute elements of Matrix U
//    // treat the first 2 U_bar and first 3 U specially
//    U[0] = LOmicron_diag[0].matrixLLT().triangularView<Lower> ().solve(-Bm_transp);
//    U_bar[0] = LPi_diag[0].matrixLLT().triangularView<Lower> ().solve(eye);
//
//    U[1] = LOmicron_diag[0].matrixLLT().triangularView<Lower> ().solve(-B_transp);
//    U_bar[1] = LRho_diag[0].matrixLLT().triangularView<Lower> ().solve(eye);
//    U[2] = LOmicron_diag[1].matrixLLT().triangularView<Lower> ().solve( /*zero*/-LSigma_offDiag[0] * U_bar[1]);
//
//    // remaining U_bar and U in middle have structures
//    for (int i = 1; i <= _N - 2; i++)
//    {
//      U_bar[2 + (i - 1) * 5] = LPi_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(-Am_tilde_transp);
//      U_bar[3 + (i - 1) * 5] = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(-Bm_bar_transp);
//      U[3 + (i - 1) * 3] = LOmicron_diag[i].matrixLLT().triangularView<Lower> ().solve(
//                                                                                       -Bm_transp - LSigma_offDiag[i
//                                                                                           - 1]
//                                                                                           * U_bar[3 + (i - 1) * 5]);
//      U_bar[4 + (i - 1) * 5] = LPi_diag[i].matrixLLT().triangularView<Lower> ().solve(eye);
//
//      U_bar[5 + (i - 1) * 5] = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(-A_bar_transp);
//      U[4 + (i - 1) * 3]
//          = LOmicron_diag[i].matrixLLT().triangularView<Lower> ().solve(
//                                                                        -B_transp - LSigma_offDiag[i - 1] * U_bar[5
//                                                                            + (i - 1) * 5]);
//      U_bar[6 + (i - 1) * 5] = LRho_diag[i].matrixLLT().triangularView<Lower> ().solve(eye);
//      U[5 + (i - 1) * 3] = LOmicron_diag[i + 1].matrixLLT().triangularView<Lower> ().solve(
//      /*zero*/-LSigma_offDiag[i] * U_bar[6 + (i - 1) * 5]);
//    }
//
//    U_bar[2 + (_N - 2) * 5] = LPi_diag[_N - 2].matrixLLT().triangularView<Lower> ().solve(-Am_tilde_transp);
//    U_bar[3 + (_N - 2) * 5] = LRho_diag[_N - 2].matrixLLT().triangularView<Lower> ().solve(-Bm_bar_transp);
//    U[3 + (_N - 2) * 3] = LOmicron_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
//                                                                                           -Bm_transp
//                                                                                               - LSigma_offDiag[_N - 2]
//                                                                                                   * U_bar[3 + (_N - 2)
//                                                                                                       * 5]);
//    U_bar[4 + (_N - 2) * 5] = LPi_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(eye);
//
//    U_bar[5 + (_N - 2) * 5] = LRho_diag[_N - 2].matrixLLT().triangularView<Lower> ().solve(-A_bar_transp);
//    U[4 + (_N - 2) * 3] = LOmicron_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
//                                                                                           -B_transp
//                                                                                               - LSigma_offDiag[_N - 2]
//                                                                                                   * U_bar[5 + (_N - 2)
//                                                                                                       * 5]);
//    U_bar[6 + (_N - 2) * 5] = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(eye);
//
//    if (_N == _pos_omega)
//      U[5 + (_N - 2) * 3] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
//      /*zero*/-LSigma_offDiag[_N - 1] * U_bar[6 + (_N - 2) * 5]);
//    else
//    {
//      UO[0] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
//      /*zero*/-LLambda0 * U_bar[1 + (_pos_omega - 1) * 5] - LLambda1 * U[2 + (_pos_omega - 1) * 3]);
//      UO[1] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
//      /*zero*/-LLambda0 * U_bar[3 + (_pos_omega - 1) * 5] - LLambda1 * U[3 + (_pos_omega - 1) * 3]);
//      UO[2] = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
//      /*zero*/-LLambda0 * U_bar[5 + (_pos_omega - 1) * 5] - LLambda1 * U[4 + (_pos_omega - 1) * 3]);
//    }
//
//    /*
//     cout << "U[0]" << endl << U[0] << endl << endl;
//     cout << "U_bar[0]" << endl << U_bar[0] << endl << endl;
//     cout << "U[1]" << endl << U[1] << endl << endl;
//     cout << "U_bar[1]" << endl << U_bar[1] << endl << endl;
//     cout << "U[2]" << endl << U[2] << endl << endl;
//     for (int i=1; i<= _N-2; i++)
//     {
//     cout << "U_bar[" << 2+(i-1)*5 << "]" << endl << U_bar[2+(i-1)*5] << endl << endl ;
//     cout << "U_bar[" << 3+(i-1)*5 << "]" << endl << U_bar[3+(i-1)*5] << endl << endl ;
//     cout << "U[" << 3+(i-1)*3 << "]" << endl << U[3+(i-1)*3] << endl << endl;
//     cout << "U_bar[" << 4+(i-1)*5 << "]" << endl << U_bar[4+(i-1)*5] << endl << endl;
//
//     cout << "U_bar[" << 5+(i-1)*5 << "]" << endl << U_bar[5+(i-1)*5] << endl << endl ;
//     cout << "U[" << 4+(i-1)*3 << "]" << endl << U[4+(i-1)*3] << endl << endl ;
//     cout << "U_bar[" << 6+(i-1)*5 << "]" << endl << U_bar[6+(i-1)*5] << endl << endl ;
//     cout << "U[" << 5+(i-1)*3 << "]" << endl << U[5+(i-1)*3] << endl << endl ;
//     }
//     cout << "U_bar[" << 2+(_N-2)*5 << "]" << endl << U_bar[2+(_N-2)*5] << endl << endl ;
//     cout << "U_bar[" << 3+(_N-2)*5 << "]" << endl << U_bar[3+(_N-2)*5] << endl << endl ;
//     cout << "U[" << 3+(_N-2)*3 << "]" << endl << U[3+(_N-2)*3] << endl << endl;
//     cout << "U_bar[" << 4+(_N-2)*5 << "]" << endl << U_bar[4+(_N-2)*5] << endl << endl;
//
//     cout << "U_bar[" << 5+(_N-2)*5 << "]" << endl << U_bar[5+(_N-2)*5] << endl << endl ;
//     cout << "U[" << 4+(_N-2)*3 << "]" << endl << U[4+(_N-2)*3] << endl << endl ;
//     cout << "U_bar[" << 6+(_N-2)*5 << "]" << endl << U_bar[6+(_N-2)*5] << endl << endl ;
//
//     if (_N == _pos_omega)
//     {
//     cout << "U[" << 5+(_N-2)*3 << "]" << endl << U[5+(_N-2)*3] << endl << endl ;
//     }
//     else
//     {
//     cout << "UO[0]" << endl << UO[0] << endl << endl;
//     cout << "UO[1]" << endl << UO[1] << endl << endl;
//     cout << "UO[2]" << endl << UO[2] << endl << endl;
//     }
//     */
//
//    // 2. Compute elements in Matrix X
//    // treat the first 2 X_bar and first 3 X specially
//
//    if (_pos_omega != _N)
//    {
//      XO[0] = LOmicron_diag_transp[_N].template triangularView<Upper> ().solve(UO[0]);
//      XO[1] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[1]);
//      XO[2] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[2]);
//    }
//
//    X[0] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[0]);
//    X_bar[0] = LPi_diag_transp[0].template triangularView<Upper>().solve(U_bar[0]);
//
//    if (_pos_omega == 1) // build col1
//
//          {
//            X[1] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[1]);
//            X[2] = LOmicron_diag_transp[1].template triangularView<Upper>().solve(U[2] - LLambda1_transp*XO[0]);
//            X_bar[1] = LRho_diag_transp[0].template triangularView<Upper>().solve(U_bar[1] - LSigma_offDiag_transp[0]*X[2] - LLambda0_transp*XO[0]);
//          }
//          else
//          {
//            X[1] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[1]);
//            X[2] = LOmicron_diag_transp[1].template triangularView<Upper>().solve(U[2]);
//            X_bar[1] = LRho_diag_transp[0].template triangularView<Upper>().solve(U_bar[1] - LSigma_offDiag_transp[0]*X[2]);
//          }
//
//          // remaining X_bar and X have structures
//          // if(_N != _pos_omega), then the off-diagonal element incluences three columns, i.e col1, col2, col3
//          for (int i = 1; i <= _N-2; i++)
//          {
//            if (i == _pos_omega) // col2
//
//            {
//              X_bar[2+(_pos_omega-1)*5] = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[2+(_pos_omega-1)*5]);
//              X[3+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(U[3+(_pos_omega-1)*3]-LLambda1_transp*XO[1]);
//              X_bar[3+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[3+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[3+(_pos_omega-1)*3] - LLambda0_transp*XO[1]);
//              X_bar[4+(_pos_omega-1)*5] = LPi_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar[4+(_pos_omega-1)*5]);
//            }
//            else // standard
//
//            {
//              X_bar[2+(i-1)*5] = LPi_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[2+(i-1)*5]);
//              X[3+(i-1)*3] = LOmicron_diag_transp[i].template triangularView<Upper>().solve(U[3+(i-1)*3]);
//              X_bar[3+(i-1)*5] = LRho_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[3+(i-1)*5] - LSigma_offDiag_transp[i-1]*X[3+(i-1)*3]);
//              X_bar[4+(i-1)*5] = LPi_diag_transp[i].template triangularView<Upper>().solve(U_bar[4+(i-1)*5]);
//            }
//
//            if (i == _pos_omega) // col3
//
//            {
//              X[4+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve( U[4+(_pos_omega-1)*3] - LLambda1_transp*XO[2]);
//              X_bar[5+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[5+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[4+(_pos_omega-1)*3] - LLambda0_transp*XO[2]);
//              X[5+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega+1].template triangularView<Upper>().solve(U[5+(_pos_omega-1)*3]);
//              X_bar[6+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar[6+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega]*X[5+(_pos_omega-1)*3]);
//            }
//            else if(i == _pos_omega-1) // col1
//
//            {
//              X[4+(_pos_omega-2)*3] = LOmicron_diag_transp[_pos_omega-1].template triangularView<Upper>().solve( U[4+(_pos_omega-2)*3] );
//              X_bar[5+(_pos_omega-2)*5] = LRho_diag_transp[_pos_omega-2].template triangularView<Upper>().solve(U_bar[5+(_pos_omega-2)*5] - LSigma_offDiag_transp[_pos_omega-2]*X[4+(_pos_omega-2)*3]);
//              X[5+(_pos_omega-2)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(U[5+(_pos_omega-2)*3] - LLambda1_transp*XO[0]);
//              X_bar[6+(_pos_omega-2)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[6+(_pos_omega-2)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[5+(_pos_omega-2)*3] - LLambda0_transp*XO[0]);
//            }
//            else // if the off-diag element in P*z<h has no influence
//
//            {
//              X[4+(i-1)*3] = LOmicron_diag_transp[i].template triangularView<Upper>().solve( U[4+(i-1)*3] );
//              X_bar[5+(i-1)*5] = LRho_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[5+(i-1)*5] - LSigma_offDiag_transp[i-1]*X[4+(i-1)*3]);
//              X[5+(i-1)*3] = LOmicron_diag_transp[i+1].template triangularView<Upper>().solve(U[5+(i-1)*3]);
//              X_bar[6+(i-1)*5] = LRho_diag_transp[i].template triangularView<Upper>().solve(U_bar[6+(i-1)*5] - LSigma_offDiag_transp[i]*X[5+(i-1)*3]);
//            }
//          }
//
//          // compute last two columns
//          if (_pos_omega == _N)
//          {
//            X_bar[2+(_N-2)*5] = LPi_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[2+(_N-2)*5]);
//            X[3+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve(U[3+(_N-2)*3]);
//            X_bar[3+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[3+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[3+(_N-2)*3]);
//            X_bar[4+(_N-2)*5] = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[4+(_N-2)*5]);
//
//            X[4+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve( U[4+(_N-2)*3] );
//            X_bar[5+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[5+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[4+(_N-2)*3]);
//            X[5+(_N-2)*3] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(U[5+(_N-2)*3]);
//            X_bar[6+(_N-2)*5] = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[6+(_N-2)*5] - LSigma_offDiag_transp[_N-1]*X[5+(_N-2)*3]);
//          }
//          else if(_pos_omega == _N-1) // compute col2 and col3
//
//          {
//            X_bar[2+(_N-2)*5] = LPi_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[2+(_N-2)*5]);
//            X[3+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve(U[3+(_N-2)*3] - LLambda1_transp*XO[1]);
//            X_bar[3+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[3+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[3+(_N-2)*3] - LLambda0_transp*XO[1]);
//            X_bar[4+(_N-2)*5] = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[4+(_N-2)*5]);
//
//            X[4+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve( U[4+(_N-2)*3] - LLambda1_transp*XO[2]);
//            X_bar[5+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[5+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[4+(_N-2)*3] - LLambda0_transp*XO[2]);
//            X_bar[6+(_N-2)*5] = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[6+(_N-2)*5]);
//          }
//          else // standard, no influence by off-diag element in P
//
//          {
//            X_bar[2+(_N-2)*5] = LPi_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[2+(_N-2)*5]);
//            X[3+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve(U[3+(_N-2)*3]);
//            X_bar[3+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[3+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[3+(_N-2)*3]);
//            X_bar[4+(_N-2)*5] = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[4+(_N-2)*5]);
//
//            X[4+(_N-2)*3] = LOmicron_diag_transp[_N-1].template triangularView<Upper>().solve( U[4+(_N-2)*3]);
//            X_bar[5+(_N-2)*5] = LRho_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar[5+(_N-2)*5] - LSigma_offDiag_transp[_N-2]*X[4+(_N-2)*3]);
//            X_bar[6+(_N-2)*5] = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar[6+(_N-2)*5]);
//          }
//
//          /*
//           cout << "X[0]" << endl << X[0] << endl << endl;
//           cout << "X_bar[0]" << endl << X_bar[0] << endl << endl;
//           cout << "X[1]" << endl << X[1] << endl << endl;
//           cout << "X_bar[1]" << endl << X_bar[1] << endl << endl;
//           cout << "X[2]" << endl << X[2] << endl << endl;
//           for (int i=1; i<= _N-2; i++)
//           {
//           cout << "X_bar[" << 2+(i-1)*5 << "]" << endl << X_bar[2+(i-1)*5] << endl << endl ;
//           cout << "X_bar[" << 3+(i-1)*5 << "]" << endl << X_bar[3+(i-1)*5] << endl << endl ;
//           cout << "X[" << 3+(i-1)*3 << "]" << endl << X[3+(i-1)*3] << endl << endl;
//           cout << "X_bar[" << 4+(i-1)*5 << "]" << endl << X_bar[4+(i-1)*5] << endl << endl;
//
//           cout << "X_bar[" << 5+(i-1)*5 << "]" << endl << X_bar[5+(i-1)*5] << endl << endl ;
//           cout << "X[" << 4+(i-1)*3 << "]" << endl << X[4+(i-1)*3] << endl << endl ;
//           cout << "X_bar[" << 6+(i-1)*5 << "]" << endl << X_bar[6+(i-1)*5] << endl << endl ;
//           cout << "X[" << 5+(i-1)*3 << "]" << endl << X[5+(i-1)*3] << endl << endl ;
//           }
//           cout << "X_bar[2+(_N-2)*5]" << endl << X_bar[2+(_N-2)*5] << endl << endl ;
//           cout << "X_bar[3+(_N-2)*5]" << endl << X_bar[3+(_N-2)*5] << endl << endl ;
//           cout << "X[3+(_N-2)*3]" << endl << X[3+(_N-2)*3] << endl << endl;
//           cout << "X_bar[4+(_N-2)*5]" << endl << X_bar[4+(_N-2)*5] << endl << endl;
//
//           cout << "X_bar[5+(_N-2)*5]" << endl << X_bar[5+(_N-2)*5] << endl << endl ;
//           cout << "X[4+(_N-2)*3]" << endl << X[4+(_N-2)*3] << endl << endl ;
//           cout << "X_bar[6+(_N-2)*5]" << endl << X_bar[6+(_N-2)*5] << endl << endl ;
//           if (_pos_omega == _N)
//           {
//           cout << "X[5+(_N-2)*3]" << endl << X[5+(_N-2)*3] << endl << endl ;
//           }
//           else
//           {
//           cout << "XO[0]" << endl << XO[0] << endl << endl;
//           cout << "XO[1]" << endl << XO[1] << endl << endl;
//           cout << "XO[2]" << endl << XO[2] << endl << endl;
//           }
//           */
//
//          // 3. Compute Y = C*X
//          // compute first three Y separately
//          Y[0][0] = -Bm*X[0] + X_bar[0];
//
//          Y[0][1] = -Bm*X[1];
//          Y[1][1] = -B*X[1] + X_bar[1];
//
//          // compute rest by filling Y column by column; done for 2 neighboring rows
//          // we fill Y column by column, treating the first two columns specially
//          for (int i=1; i <= _N-1; i++)
//          {
//            Y[0][2*(i-1)+2] = X_bar[2+(i-1)*5];
//            Y[1][2*(i-1)+2] = X_bar[3+(i-1)*5];
//            Y[2][2*(i-1)+2] = -Am_tilde*X_bar[2+(i-1)*5] - Bm_bar*X_bar[3+(i-1)*5] - Bm*X[3+(i-1)*3] + X_bar[4+(i-1)*5];
//
//            Y[0][2*(i-1)+3] = X_bar[5+(i-1)*5];
//            Y[1][2*(i-1)+3] = -Bm_bar*X_bar[5+(i-1)*5] - Bm*X[4+(i-1)*3];
//            Y[2][2*(i-1)+3] = -A_bar*X_bar[5+(i-1)*5] - B*X[4+(i-1)*3] + X_bar[6+(i-1)*5];
//          }
//
//          /*
//           cout << "Y[0][0]" << endl << Y[0][0] << endl << endl;
//           cout << "Y[0][1]" << endl << Y[0][1] << endl << endl;
//           cout << "Y[1][1]" << endl << Y[1][1] << endl << endl;
//           for (int i=1; i <= _N-1; i++)
//           {
//           cout << "Y[0][2*(i-1)+2]" << endl << Y[0][2*(i-1)+2] << endl << endl;
//           cout << "Y[1][2*(i-1)+2]" << endl << Y[1][2*(i-1)+2] << endl << endl;
//           cout << "Y[2][2*(i-1)+2]" << endl << Y[2][2*(i-1)+2] << endl << endl;
//
//           cout << "Y[0][2*(i-1)+3]" << endl << Y[0][2*(i-1)+3] << endl << endl;
//           cout << "Y[1][2*(i-1)+3]" << endl << Y[1][2*(i-1)+3] << endl << endl;
//           cout << "Y[2][2*(i-1)+3]" << endl << Y[2][2*(i-1)+3] << endl << endl;
//           }
//           */
//        }

// beta = r_C + C*(Phi\r_d_bar), where r_d_bar = -r_H - P'*(lambda.*r_P./t - r_T./t);
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compBeta()
  {
    // solve in steps
    // 1. r_d_bar = -r_H - P'*(lambda.*r_P./t - r_T./t);
    // 2. L*tmp2 = r_d_bar   --> compute tmp2
    // 3. L'*tmp1 = tmp2  ---> compute tmp1
    // 4. beta = r_C + C*tmp1, tmp1 = Phi\r_d_bar

    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> tmp; // tmp = lambda.*r_P./t - r_T./t

    // 1. r_d_bar = -r_H - P'*(lambda.*r_P./t - r_T./t)
    tmp = (lambda.cwiseProduct(r_P) - r_T).cwiseQuotient(slack);
    // 1.1 r_d_bar = -r_H - P'*tmp
    r_d_bar.template segment<_m> (0) = -r_H.template segment<_m> (0) - Fu_transp[0] * tmp.template segment<_nInp> (0);
    int offset1 = 2 * _n; // offset required in nu for C'*nu
    for (int i = 1; i <= _N - 1; i++)
    {
      r_d_bar.template segment<_n> (_m + (i - 1) * offset) = -r_H.template segment<_n> (_m + (i - 1) * offset);

      if (i != _pos_omega)
      {
        r_d_bar.template segment<_n> (_m + (i - 1) * offset + _n) = -r_H.template segment<_n> (
                                                                                               _m + (i - 1) * offset
                                                                                                   + _n) - Fx_transp[i
            - 1] * tmp.template segment<_nSt> (_nInp + (i - 1) * (_nInp + _nSt)) - Fu_bar_transp[i]
            * tmp.template segment<_nInp> (i * (_nInp + _nSt));
      }
      else // must add the additional term: F_xTheta'*d(.)
      {
        r_d_bar.template segment<_n> (_m + (_pos_omega - 1) * offset + _n) = -r_H.template segment<_n> (
                                                                                                        _m + (i - 1)
                                                                                                            * offset
                                                                                                            + _n)
            - Fx_transp[_pos_omega - 1] * tmp.template segment<_nSt> (_nInp + (i - 1) * (_nInp + _nSt))
            - Fu_bar_transp[_pos_omega] * tmp.template segment<_nInp> (i * (_nInp + _nSt)) - F_xTheta_transp
            * tmp.template segment<_nF_xTheta> (_N * (_nSt + _nInp)); // used to be: num_constr-_nF_xTheta;
      }

      r_d_bar.template segment<_m> (_m + (i - 1) * offset + _n + _n) = -r_H.template segment<_m> (
                                                                                                  _m + (i - 1) * offset
                                                                                                      + _n + _n)
          - Fu_transp[i] * tmp.template segment<_nInp> (i * (_nInp + _nSt));

    }
    r_d_bar.template segment<_n> (_m + (_N - 1) * offset) = -r_H.template segment<_n> (_m + (_N - 1) * offset);

    if (_pos_omega == _N)
    {
      r_d_bar.template segment<_n> (_m + (_N - 1) * offset + _n) = -r_H.template segment<_n> (
                                                                                              _m + (_N - 1) * offset
                                                                                                  + _n) - Fx_transp[_N
          - 1] * tmp.template segment<_nSt> (_nInp + (_N - 1) * (_nInp + _nSt)) - F_xTheta_transp
          * tmp.template segment<_nF_xTheta> (_N * (_nInp + _nSt));
    }
    else //standard
    {
      r_d_bar.template segment<_n> (_m + (_N - 1) * offset + _n) = -r_H.template segment<_n> (
                                                                                              _m + (_N - 1) * offset
                                                                                                  + _n) - Fx_transp[_N
          - 1] * tmp.template segment<_nSt> (_nInp + (_N - 1) * (_nInp + _nSt));
    }
    r_d_bar.template segment<_m> (_m + (_N - 1) * offset + _n + _n) = -r_H.template segment<_m> (
                                                                                                 _m + (_N - 1) * offset
                                                                                                     + _n + _n)
        - F_theta_transp * tmp.template segment<_nF_xTheta> (_N * (_nSt + _nInp));
    // cout << setprecision(30) << "r_d_bar" << endl << r_d_bar << endl << endl;


    // 2. compute tmp2: L*tmp2 = r_d_bar
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> tmp2;
    tmp2.template segment<_m> (0)
        = LOmicron_diag[0].matrixLLT().triangularView<Lower> ().solve(r_d_bar.template segment<_m> (0));
    for (int i = 1; i <= _N - 1; i++)
    {
      tmp2.template segment<_n> (_m + (i - 1) * offset)
          = LPi_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                       r_d_bar.template segment<_n> (
                                                                                                     _m + (i - 1)
                                                                                                         * offset));
      tmp2.template segment<_n> (_m + (i - 1) * offset + _n)
          = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                        r_d_bar.template segment<_n> (
                                                                                                      _m + (i - 1)
                                                                                                          * offset + _n));
      tmp2.template segment<_m> (_m + (i - 1) * offset + _n + _n)
          = LOmicron_diag[i].matrixLLT().triangularView<Lower> ().solve(
                                                                        r_d_bar.template segment<_m> (
                                                                                                      _m + (i - 1)
                                                                                                          * offset + _n
                                                                                                          + _n)
                                                                            - LSigma_offDiag[i - 1]
                                                                                * tmp2.template segment<_n> (
                                                                                                             _m
                                                                                                                 + (i
                                                                                                                     - 1)
                                                                                                                     * offset
                                                                                                                 + _n));
    }
    if (_pos_omega == _N)
    {
      tmp2.template segment<_n> (_m + (_N - 1) * offset)
          = LPi_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                        r_d_bar.template segment<_n> (
                                                                                                      _m + (_N - 1)
                                                                                                          * offset));
      tmp2.template segment<_n> (_m + (_N - 1) * offset + _n)
          = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                         r_d_bar.template segment<_n> (
                                                                                                       _m + (_N - 1)
                                                                                                           * offset
                                                                                                           + _n));
      tmp2.template segment<_m> (_m + (_N - 1) * offset + _n + _n)
          = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
                                                                         r_d_bar.template segment<_m> (
                                                                                                       _m + (_N - 1)
                                                                                                           * offset
                                                                                                           + _n + _n)
                                                                             - LSigma_offDiag[_N - 1]
                                                                                 * tmp2.template segment<_n> (
                                                                                                              _m
                                                                                                                  + (_N
                                                                                                                      - 1)
                                                                                                                      * offset
                                                                                                                  + _n));
    }
    else // need to use LLambda0 and LLambda1
    {
      tmp2.template segment<_n> (_m + (_N - 1) * offset)
          = LPi_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                        r_d_bar.template segment<_n> (
                                                                                                      _m + (_N - 1)
                                                                                                          * offset));
      tmp2.template segment<_n> (_m + (_N - 1) * offset + _n)
          = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                         r_d_bar.template segment<_n> (
                                                                                                       _m + (_N - 1)
                                                                                                           * offset
                                                                                                           + _n));
      tmp2.template segment<_m> (_m + (_N - 1) * offset + _n + _n)
          = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
                                                                         r_d_bar.template segment<_m> (
                                                                                                       _m + (_N - 1)
                                                                                                           * offset
                                                                                                           + _n + _n)
                                                                             - LLambda0
                                                                                 * tmp2.template segment<_n> (
                                                                                                              _m
                                                                                                                  + (_pos_omega
                                                                                                                      - 1)
                                                                                                                      * offset
                                                                                                                  + _n)
                                                                             - LLambda1
                                                                                 * tmp2.template segment<_m> (
                                                                                                              _m
                                                                                                                  + (_pos_omega
                                                                                                                      - 1)
                                                                                                                      * offset
                                                                                                                  + _n
                                                                                                                  + _n));
    }
    // cout << "tmp2:" << setprecision(30)<< endl << tmp2 << endl << endl;


    // 3. compute tmp1: L'*tmp1 = tmp2
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> tmp1;
    tmp1.template segment<_m> (0)
        = LOmicron_diag_transp[0].template triangularView<Upper> ().solve(tmp2.template segment<_m>(0));
        for (int i = 1; i <= _pos_omega-1; i++)
        {
          tmp1.template segment<_n>(_m+(i-1)*offset) = LPi_diag_transp[i-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(i-1)*offset));
          tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n) = LOmicron_diag_transp[i].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(i-1)*offset+_n+_n));
          tmp1.template segment<_n>(_m+(i-1)*offset+_n) = LRho_diag_transp[i-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(i-1)*offset+_n) - LSigma_offDiag_transp[i-1]*tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n));
        }
        // the missing block is computed after the last block is computed
        for (int i = _pos_omega+1; i <= _N-1; i++)
        {
          tmp1.template segment<_n>(_m+(i-1)*offset) = LPi_diag_transp[i-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(i-1)*offset));
          tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n) = LOmicron_diag_transp[i].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(i-1)*offset+_n+_n));
          tmp1.template segment<_n>(_m+(i-1)*offset+_n) = LRho_diag_transp[i-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(i-1)*offset+_n) - LSigma_offDiag_transp[i-1]*tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n));
        }
        // last block
        if (_pos_omega == _N)
        {
          tmp1.template segment<_n>(_m+(_N-1)*offset) = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_N-1)*offset));
          tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(_N-1)*offset+_n+_n));
          tmp1.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_N-1)*offset+_n) - LSigma_offDiag_transp[_N-1]*tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n));
        }
        else // standard and compute missing block

        {
          tmp1.template segment<_n>(_m+(_N-1)*offset) = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_N-1)*offset));
          tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(_N-1)*offset+_n+_n));
          tmp1.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_N-1)*offset+_n));

          tmp1.template segment<_n>(_m+(_pos_omega-1)*offset) = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_pos_omega-1)*offset));
          tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda1_transp*tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) );
          tmp1.template segment<_n>(_m+(_pos_omega-1)*offset+_n) = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_pos_omega-1)*offset+_n) - LSigma_offDiag_transp[_pos_omega-1]*tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda0_transp*tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) );
        }
        // cout << setprecision(30) << "tmp1:" << endl << tmp1 << endl << endl;

        // cout << setprecision(30) << "r_C" << endl << r_C << endl << endl;
        // 4. beta = r_C + C*tmp1
        beta.template segment<_n>(0) = r_C.template segment<_n>(0) + ( -Bm*tmp1.template segment<_m>(0) + tmp1.template segment<_n>(_m) );
        beta.template segment<_n>(_n) = r_C.template segment<_n>(_n) + (- B*tmp1.template segment<_m>(0) + tmp1.template segment<_n>(_m+_n) );

        for(int i=1; i<= _N-1; i++)
        {
          beta.template segment<_n>(2*i*_n) = r_C.template segment<_n>(2*i*_n) + (
              - Am_tilde*tmp1.template segment<_n>(_m+(i-1)*offset)
              - Bm_bar * tmp1.template segment<_n>(_m+(i-1)*offset+_n)
              - Bm* tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n)
              + tmp1.template segment<_n>(i*offset+_m) );

          beta.template segment<_n>(2*i*_n+_n) = r_C.template segment<_n>( 2*i*_n+_n) + (
              - A_bar * tmp1.template segment<_n>(_m+_n+(i-1)*offset)
              - B * tmp1.template segment<_m>(_m+_n+(i-1)*offset+_n)
              + tmp1.template segment<_n>(_m+(i)*offset+_n) );
        }
        // cout << setprecision(30) << "beta:" << endl << beta << endl << endl;
      }

    // ------------ function computes L: L*L' = Y --------------------------
    // -------------- compute components of L[3][2*_N] ---------------------
    // ---- remark: L[i][i] are lower triangular matrices and can be accessed by L[i][i].matrixLLT().triangularView<Lower>()
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compL() // L*L' = Y
  {
    // special treatment for L11, L21, L31
    L_diag[0].compute(Y[0][0]); // Y11 = L11*L11'
    L_diag_transp[0] = L_diag[0].matrixLLT().transpose();

    L_offDiag_transp[0][0] = L_diag[0].matrixLLT().triangularView<Lower> ().solve(Y[0][1]);
    L_offDiag[0][0] = L_offDiag_transp[0][0].transpose();

    L_offDiag_transp[1][0] = L_diag[0].matrixLLT().triangularView<Lower> ().solve(Y[0][2]);
    L_offDiag[1][0] = L_offDiag_transp[1][0].transpose();

    // special treatment for L22, L32, L42
    L_diag[1].compute(Y[1][1] - L_offDiag[0][0] * L_offDiag_transp[0][0]);
    L_diag_transp[1] = L_diag[1].matrixLLT().transpose();

    L_offDiag_transp[0][1] = L_diag[1].matrixLLT().triangularView<Lower> ().solve(
                                                                                  Y[1][2] - L_offDiag[0][0]
                                                                                      * L_offDiag_transp[1][0]);
    L_offDiag[0][1] = L_offDiag_transp[0][1].transpose();

    L_offDiag_transp[1][1] = L_diag[1].matrixLLT().triangularView<Lower> ().solve(Y[0][3]);
    L_offDiag[1][1] = L_offDiag_transp[1][1].transpose();

    // cases in the middle
    for (int i = 1; i <= 2 * _N - 4; i++)
    {
      L_diag[i + 1].compute(
                            Y[2][i + 1] - L_offDiag[1][i - 1] * L_offDiag_transp[1][i - 1] - L_offDiag[0][i]
                                * L_offDiag_transp[0][i]);
      L_diag_transp[i + 1] = L_diag[i + 1].matrixLLT().transpose();
      L_offDiag_transp[0][i + 1]
          = L_diag[i + 1].matrixLLT().triangularView<Lower> ().solve(
                                                                     Y[1][i + 2] - L_offDiag[0][i]
                                                                         * L_offDiag_transp[1][i]);
      L_offDiag[0][i + 1] = L_offDiag_transp[0][i + 1].transpose();

      L_offDiag_transp[1][i + 1] = L_diag[i + 1].matrixLLT().triangularView<Lower> ().solve(Y[0][i + 3]);
      L_offDiag[1][i + 1] = L_offDiag_transp[1][i + 1].transpose();
    }

    // special treatment in the end, i.e. i = 2*_N-3
    L_diag[2 * _N - 2].compute(
                               Y[2][2 * _N - 2] - L_offDiag[1][2 * _N - 4] * L_offDiag_transp[1][2 * _N - 4]
                                   - L_offDiag[0][2 * _N - 3] * L_offDiag_transp[0][2 * _N - 3]);
    L_diag_transp[2 * _N - 2] = L_diag[2 * _N - 2].matrixLLT().transpose();
    L_offDiag_transp[0][2 * _N - 2]
        = L_diag[2 * _N - 2].matrixLLT().triangularView<Lower> ().solve(
                                                                        Y[1][2 * _N - 1] - L_offDiag[0][2 * _N - 3]
                                                                            * L_offDiag_transp[1][2 * _N - 3]);
    L_offDiag[0][2 * _N - 2] = L_offDiag_transp[0][2 * _N - 2].transpose();

    // i = 2*_N-2
    L_diag[2 * _N - 1].compute(
                               Y[2][2 * _N - 1] - L_offDiag[1][2 * _N - 3] * L_offDiag_transp[1][2 * _N - 3]
                                   - L_offDiag[0][2 * _N - 2] * L_offDiag_transp[0][2 * _N - 2]);
    L_diag_transp[2 * _N - 1] = L_diag[2 * _N - 1].matrixLLT().transpose();
  }

// ------------ function computes L*L'*dnu = -beta ---------------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compDnu()
  {
    // 1) first, solve for delta: L*delta = -beta
    Matrix<Type, 2 * _N * _n, 1> delta;

    // special cases in the beginning
    delta.template segment<_n> (0)
        = L_diag[0].matrixLLT().triangularView<Lower> ().solve(-beta.template segment<_n> (0));
    delta.template segment<_n> (_n)
        = L_diag[1].matrixLLT().triangularView<Lower> ().solve(
                                                               -beta.template segment<_n> (_n) - L_offDiag[0][0]
                                                                   * delta.template segment<_n> (0));

    // remaining cases are regular
    for (int i = 1; i <= 2 * _N - 2; i++)
      delta.template segment<_n> (_n + i * _n)
          = L_diag[i + 1].matrixLLT().triangularView<Lower> ().solve(
                                                                     -beta.template segment<_n> (_n + i * _n)
                                                                         - L_offDiag[1][i - 1]
                                                                             * delta.template segment<_n> ((i - 1) * _n)
                                                                         - L_offDiag[0][i]
                                                                             * delta.template segment<_n> (i * _n));

    // 2) now, solve for L'*Dnu = delta
    dnu.template segment<_n> (2 * _n * _N - _n)
        = L_diag_transp[2 * _N - 1].template triangularView<Upper> ().solve(delta.template segment<_n>(2*_n*_N - _n) );
        dnu.template segment<_n>(2*_n*_N - _n - _n) = L_diag_transp[2*_N-2].template triangularView<Upper>().solve( delta.template segment<_n>(2*_n*_N - _n - _n) - L_offDiag_transp[0][2*_N-2]*dnu.template segment<_n>(2*_n*_N - _n) );

        //remaining cases are regular
        for (int i=1; i<=2*_N-2; i++)
        dnu.template segment<_n>(2*_n*_N-(i+2)*_n) = L_diag_transp[2*_N-(i+2)].template triangularView<Upper>().solve( delta.template segment<_n>(2*_n*_N-(i+2)*_n) - L_offDiag_transp[0][2*_N-(i+2)]*dnu.template segment<_n>(2*_n*_N-(i+1)*_n) - L_offDiag_transp[1][2*_N-(i+2)]*dnu.template segment<_n>(2*_n*_N-i*_n) );
        //cout << "dnu" << endl << dnu << endl << endl;
      }

    // ------------ function computes Phi*dz = -r_d_bar - C'*dnu ---------------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compDz()
  {
    // computed in two parts
    // 1. tmp = -r_d - C' * dnu
    // 2. L_Phi*L_Phi'*dz = tmp
    // 3. L_Phi*tmp1 = tmp
    // 4. L_Phi'*dz = tmp1

    Matrix<Type, _N * (_m + _n + _n) + _m, 1> tmp;
    tmp.template segment<_m> (0) = -r_d_bar.template segment<_m> (0) + Bm_transp * dnu.template segment<_n> (0)
        + B_transp * dnu.template segment<_n> (_n);

    for (int i = 1; i <= _N - 1; i++)
    {
      tmp.template segment<_n> (_m + (i - 1) * offset) = -r_d_bar.template segment<_n> (_m + (i - 1) * offset)
          - dnu.template segment<_n> (2 * (i - 1) * _n) + Am_tilde_transp * dnu.template segment<_n> (
                                                                                                      2 * (i - 1) * _n
                                                                                                          + _n + _n);
      tmp.template segment<_n> (_m + (i - 1) * offset + _n)
          = -r_d_bar.template segment<_n> (_m + (i - 1) * offset + _n) - dnu.template segment<_n> (
                                                                                                   2 * (i - 1) * _n
                                                                                                       + _n)
              + Bm_bar_transp * dnu.template segment<_n> (2 * (i - 1) * _n + _n + _n) + A_bar_transp
              * dnu.template segment<_n> (2 * (i - 1) * _n + _n + _n + _n);
      tmp.template segment<_m> (_m + (i - 1) * offset + _n + _n) = -r_d_bar.template segment<_m> (
                                                                                                  _m + (i - 1) * offset
                                                                                                      + _n + _n)
          + Bm_transp * dnu.template segment<_n> (2 * (i - 1) * _n + _n + _n) + B_transp
          * dnu.template segment<_n> (2 * (i - 1) * _n + _n + _n + _n);
    }

    tmp.template segment<_n> (_m + (_N - 1) * offset) = -r_d_bar.template segment<_n> (_m + (_N - 1) * offset)
        - dnu.template segment<_n> (2 * (_N - 1) * _n);
    tmp.template segment<_n> (_m + (_N - 1) * offset + _n)
        = -r_d_bar.template segment<_n> (_m + (_N - 1) * offset + _n) - dnu.template segment<_n> (
                                                                                                  2 * (_N - 1) * _n
                                                                                                      + _n);
    tmp.template segment<_m> (_m + (_N - 1) * offset + _n + _n) = -r_d_bar.template segment<_m> (
                                                                                                 _m + (_N - 1) * offset
                                                                                                     + _n + _n);
    //cout << "tmp:" << endl << tmp << endl << endl;


    // 3. L_Phi*tmp1 = tmp
    Matrix<Type, _N * (_m + _n + _n) + _m, 1> tmp1;
    tmp1.template segment<_m> (0)
        = LOmicron_diag[0].matrixLLT().triangularView<Lower> ().solve(tmp.template segment<_m> (0));
    for (int i = 1; i <= _N - 1; i++)
    {
      tmp1.template segment<_n> (_m + (i - 1) * offset)
          = LPi_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(tmp.template segment<_n> (_m + (i - 1) * offset));
      tmp1.template segment<_n> (_m + (i - 1) * offset + _n)
          = LRho_diag[i - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                        tmp.template segment<_n> (
                                                                                                  _m + (i - 1) * offset
                                                                                                      + _n));
      tmp1.template segment<_m> (_m + (i - 1) * offset + _n + _n)
          = LOmicron_diag[i].matrixLLT().triangularView<Lower> ().solve(
                                                                        tmp.template segment<_m> (
                                                                                                  _m + (i - 1) * offset
                                                                                                      + _n + _n)
                                                                            - LSigma_offDiag[i - 1]
                                                                                * tmp1.template segment<_n> (
                                                                                                             _m
                                                                                                                 + (i
                                                                                                                     - 1)
                                                                                                                     * offset
                                                                                                                 + _n));
    }
    if (_pos_omega == _N)
    {
      tmp1.template segment<_n> (_m + (_N - 1) * offset)
          = LPi_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                        tmp.template segment<_n> (
                                                                                                  _m + (_N - 1)
                                                                                                      * offset));
      tmp1.template segment<_n> (_m + (_N - 1) * offset + _n)
          = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                         tmp.template segment<_n> (
                                                                                                   _m + (_N - 1)
                                                                                                       * offset + _n));
      tmp1.template segment<_m> (_m + (_N - 1) * offset + _n + _n)
          = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
                                                                         tmp.template segment<_m> (
                                                                                                   _m + (_N - 1)
                                                                                                       * offset + _n
                                                                                                       + _n)
                                                                             - LSigma_offDiag[_N - 1]
                                                                                 * tmp1.template segment<_n> (
                                                                                                              _m
                                                                                                                  + (_N
                                                                                                                      - 1)
                                                                                                                      * offset
                                                                                                                  + _n));
    }
    else
    {
      tmp1.template segment<_n> (_m + (_N - 1) * offset)
          = LPi_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                        tmp.template segment<_n> (
                                                                                                  _m + (_N - 1)
                                                                                                      * offset));
      tmp1.template segment<_n> (_m + (_N - 1) * offset + _n)
          = LRho_diag[_N - 1].matrixLLT().triangularView<Lower> ().solve(
                                                                         tmp.template segment<_n> (
                                                                                                   _m + (_N - 1)
                                                                                                       * offset + _n));
      tmp1.template segment<_m> (_m + (_N - 1) * offset + _n + _n)
          = LOmicron_diag[_N].matrixLLT().triangularView<Lower> ().solve(
                                                                         tmp.template segment<_m> (
                                                                                                   _m + (_N - 1)
                                                                                                       * offset + _n
                                                                                                       + _n)
                                                                             - LLambda0
                                                                                 * tmp1.template segment<_n> (
                                                                                                              _m
                                                                                                                  + (_pos_omega
                                                                                                                      - 1)
                                                                                                                      * offset
                                                                                                                  + _n)
                                                                             - LLambda1
                                                                                 * tmp1.template segment<_m> (
                                                                                                              _m
                                                                                                                  + (_pos_omega
                                                                                                                      - 1)
                                                                                                                      * offset
                                                                                                                  + _n
                                                                                                                  + _n));
    }
    // cout << "tmp1:" << endl << tmp1 << endl << endl;


    // 4. L_Phi'*dz = tmp1
    dz.template segment<_m> (0)
        = LOmicron_diag_transp[0].template triangularView<Upper> ().solve(tmp1.template segment<_m>(0));
        for (int i = 1; i <= _pos_omega-1; i++)
        {
          dz.template segment<_n>(_m+(i-1)*offset) = LPi_diag_transp[i-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(i-1)*offset));
          dz.template segment<_m>(_m+(i-1)*offset+_n+_n) = LOmicron_diag_transp[i].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n));
          dz.template segment<_n>(_m+(i-1)*offset+_n) = LRho_diag_transp[i-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(i-1)*offset+_n) - LSigma_offDiag_transp[i-1]*dz.template segment<_m>(_m+(i-1)*offset+_n+_n));
        }
        // missing block is computed together with last block
        for (int i = _pos_omega+1; i <= _N-1; i++)
        {
          dz.template segment<_n>(_m+(i-1)*offset) = LPi_diag_transp[i-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(i-1)*offset));
          dz.template segment<_m>(_m+(i-1)*offset+_n+_n) = LOmicron_diag_transp[i].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n));
          dz.template segment<_n>(_m+(i-1)*offset+_n) = LRho_diag_transp[i-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(i-1)*offset+_n) - LSigma_offDiag_transp[i-1]*dz.template segment<_m>(_m+(i-1)*offset+_n+_n));
        }

        // last block
        if (_pos_omega == _N)
        {
          dz.template segment<_n>(_m+(_N-1)*offset) = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_N-1)*offset));
          dz.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n));
          dz.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_N-1)*offset+_n) - LSigma_offDiag_transp[_N-1]*dz.template segment<_m>(_m+(_N-1)*offset+_n+_n));
        }
        else // standard ending and missing block

        {
          dz.template segment<_n>(_m+(_N-1)*offset) = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_N-1)*offset));
          dz.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n));
          dz.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_N-1)*offset+_n) );

          dz.template segment<_n>(_m+(_pos_omega-1)*offset) = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_pos_omega-1)*offset));
          dz.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda1_transp*dz.template segment<_m>(_m+(_N-1)*offset+_n+_n) );
          dz.template segment<_n>(_m+(_pos_omega-1)*offset+_n) = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_pos_omega-1)*offset+_n) - LSigma_offDiag_transp[_pos_omega-1]*dz.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda0_transp*dz.template segment<_m>(_m+(_N-1)*offset+_n+_n) );
        }
        // cout << setprecision(30) << "dz" << endl << dz << endl << endl;
      }

    // ------------ function computes dlambda = (P*dz - r_P + r_T./lambda) .* lambda ./ t ---------------------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compDlambda()
  {
    Matrix<Type, _N * (_nSt + _nInp) + _nF_xTheta, 1> tmp; // tmp = P*dz
    // special treatment at beginning
    c_tmp = dz.template segment<_m> (0);
    tmp.template segment<_nInp> (0) = (Fu[0] * c_tmp); // should be >0
    // general treatment in the middle, class variable: offset = _n + _n + _m
    for (int i = 1; i <= _N - 1; i++) // blocks in the middle
    { // compute (h - P*z)_i
      x_bar_tmp = dz.template segment<_n> ((i - 1) * offset + _m + _n);
      c_tmp = dz.template segment<_m> ((i - 1) * offset + _m + _n + _n);
      tmp.template segment<_nSt> (_nInp + (i - 1) * (_nInp + _nSt)) = Fx[i - 1] * x_bar_tmp;
      tmp.template segment<_nInp> (i * (_nInp + _nSt)) = (Fu[i] * K * x_bar_tmp + Fu[i] * c_tmp);
    }
    // special case for last blocks
    x_bar_tmp = dz.template segment<_n> ((_N - 1) * offset + _m + _n);
    c_tmp = dz.template segment<_m> ((_N - 1) * offset + _m + _n + _n);
    tmp.template segment<_nSt> (_nInp + (_N - 1) * (_nSt + _nInp)) = (Fx[_N - 1] * x_bar_tmp);
    x_bar_tmp = dz.template segment<_n> ((_pos_omega - 1) * offset + _m + _n);
    tmp.template segment<_nF_xTheta> (_N * (_nSt + _nInp)) = (F_xTheta * x_bar_tmp + F_theta * c_tmp);
    // cout << setprecision(30) << "tmp" << endl << tmp << endl << endl;

    dlambda = (lambda.cwiseProduct(tmp - r_P) + r_T).cwiseQuotient(slack);
    // cout << setprecision(30) << "dlambda" << endl << dlambda << endl << endl;
  }

// ------------ function looks for step length (predictor and final) alpha --------
// need lambda and slack stay positive for the step size
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compAlpha()
  {
    alpha = 1;
    double alpha_tmp;
    for (int i = 0; i <= _N * (_nSt + _nInp) + _nF_xTheta - 1; i++)
    {
      if (dlambda[i] < 0)
      {
        alpha_tmp = -lambda[i] / dlambda[i];
        alpha = alpha_tmp < alpha ? alpha_tmp : alpha;
      }
      if (dslack[i] < 0)
      {
        alpha_tmp = -slack[i] / dslack[i];
        alpha = alpha_tmp < alpha ? alpha_tmp : alpha;
      }
    }
    // cout << "alpha: " << alpha << endl << endl;
  }

// ------------ function computes cost vectors q_tilde_vec, q_bar_vec, r_vec --------
template<class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
  void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::compRQ()
  {
    // compute the vectors u_star
    // "ColPivHouseholderQRPreconditioner" is more accurate, "HouseholderQRPreconditioner" is faster
    JacobiSVD<MatrixXd, HouseholderQRPreconditioner> svd(Bm, ComputeThinU | ComputeThinV); // computes SVD of Bm
    u_star[0] = svd.solve(x_star[0] - Am_tilde * (*x_hat) - tm_tilde);

    for (int i = 1; i <= _N - 1; i++)
      u_star[i] = svd.solve(x_star[i] - Am_tilde * x_star[i - 1] - tm_tilde);

    /*
     for (int i = 0; i <= _N-1; i++)
     cout << "u_star[" << i << "]" << endl << u_star[i] << endl << endl;
     */

    // compute the vectors q_bar_vec[]; q_tilde_vec[]; r_vec[]
    q_tilde_vec[0] = -2 * Q_tilde * x_star[0];
    r_vec[0] = -2 * R * u_star[0];
    // q_bar_vec[0] is never used
    for (int i = 1; i <= _N - 2; i++) // be careful how to use it
    {
      q_tilde_vec[i] = -2 * Q_tilde * x_star[i];
      r_vec[i] = -2 * R * u_star[i];
      q_bar_vec[i] = K_transp * r_vec[i];
    }
    q_tilde_vec[_N - 1] = -2 * Q_tilde_f * x_star[_N - 1];
    r_vec[_N - 1] = -2 * R * u_star[_N - 1];
    q_bar_vec[_N - 1] = K_transp * r_vec[_N - 1];

    /*
     for (int i = 0; i <= _N-1; i++)
     {
     cout << "q_tilde_vec[" << i << "]" << endl << q_tilde_vec[i] << endl << endl;
     cout << "r_vec[" << i << "]" << endl << r_vec[i] << endl << endl;
     cout << "q_bar_vec[" << i << "]" << endl << q_bar_vec[i] << endl << endl;
     }
     */
  }

#endif
