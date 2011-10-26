// Ompctp.h -- a Learning Based MPC class template
// header for class LBmpcTP
// date: September 6, 2011
// author: Xiaojing ZHANG
// version: 0.0


/* Remarks:

2) horizon >= 3 assumed
3) issues with LLT, L_diag[2*Horizon] for both Y and Phi --> save it into normal matrices
4) issues of defining length and Type for vector d, use diagonal Matrix for d^2
5) Optimize index calculations, merge different loops into 1 loop
6) avoid using namespace Eigen/std
9) upper bound size of d_diag in compPhi and compPhi_hat_PhaseI()
10) omit isFeasible in step() and phaseI() and testPos() and testPos_hat();
11) compGamma(): use vector addition rather than piece by piece addition?
12) don't need testPos() after having computed T
13) replace offset and offset1
14) move compMatrices_PhaseI() into constructor to avoid initialization overhead when calling step()
15) use DiagonalMatrix for eye, esp. PhaseI()
*/


#ifndef LBMPCTP_H_
#define LBMPCTP_H_

#include <iostream>
#include <Eigen/Dense>
#include <iomanip>		// for setprecision (9)

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

template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
class LBmpcTP
{
private:
// ----------- input variables ------------
	double kappa;	// barrier parameter
	double kappa_start;	// each new MPC step in PhaseII is initialized with kappa_start
	int n_iter;		// number of Newton iterations
	double mu;		// decrease rate of kappa, i.e. kappa := mu* kappa, mu < 1
	double eps_barrier;	// accept eps_barrier - suboptimal solution of original problem,
						// given by: kappa * num_constr < eps_barrier
	double eps_nt;	// tolerance for residua norm(r_p, r_d) in Newton method
	double eps_ntSq;
	double eps_normRp;	// tolerance for norm on primal residuum, i.e. for a fixed kappa, we want to achieve C*z = b
	double eps_normRpSq;
	double eps_ls;

	int num_constr;		// counts number of constraints for [kappa * num_constr < eps_barrier]
	double t;			// step size in backtracking line search, automatically handled
	double alpha_ls;	// 0 < alpha < 0.5 ; parameter in backtracking line search; usually 0.01 < alpha_ls < 0.3
	double beta_ls;		// 0 < alpha < 1; parameter in backtracking line search; usually 0.1 < beta_ls < 0.8
	int offset;			// = _n + _n + _m: often used variable 
	
	Matrix<Type, _n, _n> A;			// dynamics matrix A, n x n
	Matrix<Type, _n, _n> A_transp;	// transpose of above
	Matrix<Type, _n, _m> B;			// input matrix B, n x m
	Matrix<Type, _m, _n> B_transp;	// transpose of above
	Matrix<Type, _n, _n> Q_tilde;	// cost matrix on state x_tilde, p.d., n x n
	Matrix<Type, _n, _n> Q_tilde_f;	// final cost matrix, p.d.
	Matrix<Type, _m, _m> R;			// cost matrix on input, p.d. if Fu not full rank
	Matrix<Type, _n, 1> q_tilde_vec[_N];	// cost matrix
	Matrix<Type, _m, 1> r_vec[_N];	
	Matrix<Type, _nSt, _n> Fx[_N];	// array of full-rank state constraint matrices
												// #rows may change, but <= nSt
	Matrix<Type, _n, _nSt> Fx_transp[_N];	// transpose of above
	Matrix<Type, _nSt, 1> fx[_N];	// array containing right-hand-side of state constraints
												// at most of dimension nSt
	Matrix<Type, _nInp, _m> Fu[_N];	// array of (full-rank) input constraint matrices
												// #row may change, but <= nInp
	Matrix<Type, _m, _nInp> Fu_transp[_N];	//transpose of above
	Matrix<Type, _nInp, 1> fu[_N];	// array containing right-hand-side of input constraints
												// each element of at most dimension nInp
	Matrix<Type, _nF_xTheta, _n> F_xTheta;	// F_xTheta*x[m+N]+F_theta*theta <= f_xTheta
	Matrix<Type, _n, _nF_xTheta> F_xTheta_transp;	// transpose of above
	Matrix<Type, _nF_xTheta, _m> F_theta;	// full-rank constraint matrix on theta
	Matrix<Type, _m, _nF_xTheta> F_theta_transp;	// transpose of above 
	Matrix<Type, _nF_xTheta, 1> f_xTheta;	// right-hand-side of constraint above
	
	Matrix<Type, _n, _n> Lm;		// oracle matrix, O = Lm*x_tilde + Mm*u + tm
	Matrix<Type, _n, _m> Mm;		// oracle matrix, O = Lm*x_tilde + Mm*u + tm
	Matrix<Type, _n, 1> tm;			// oracle matrix, O = Lm*x_tilde + Mm*u + tm
	Matrix<Type, _m, _n> K;			// Gain matrix, s.t. (A+B*K) is stable
	Matrix<Type, _n, _m> K_transp;
	
	Matrix<Type, _n, 1> x_hat;		// state estimate of current state x_hat[m]
	Matrix<Type, _n, 1> s;			// affine offsets in linear dynamics
	Matrix<Type, _n, 1> tm_tilde;	// tm_tilde = s + tm
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> z_warm;	// warm-start z-vector, with optimal offset c[m]

	
// ----------- auxiliary variables ------------
	Matrix<Type, _n, _m> S;		// S = K'*R
	Matrix<Type, _m, _n> S_transp;	// S_trans = S'
	Matrix<Type, _n, _n> Q_bar;	// Q_bar = K'*R*K = S*K
	Matrix<Type, _n, 1> q_bar;	// q_bar' = r'*K;
	Matrix<Type, _n, 1> q_bar_vec[_N];
	Matrix<Type, _n, _n> A_bar;	// A_bar = A+B*K
	Matrix<Type, _n, _n> A_bar_transp;	// transpose of above
	Matrix<Type, _n, _n> Am_tilde;	// Am_tilde = A + Lm
	Matrix<Type, _n, _n> Am_tilde_transp;	// transpose of above
	Matrix<Type, _n, _m> Bm;	// Bm = B + Mm
	Matrix<Type, _m, _n> Bm_transp;	// transpose of above
	Matrix<Type, _n, _n> Bm_bar;	// Bm_bar = Bm * K
	Matrix<Type, _n, _n> Bm_bar_transp;	// transpose of above
	Matrix<Type, _nInp, _n> Fu_bar[_N];	// array of (full-rank) input constraint matrices, #row may change, but <= nInp
	Matrix<Type, _n, _nInp> Fu_bar_transp[_N];	//transpose of above
	
	// used in testPos(), 
	Matrix<Type, _m, 1> c_tmp;		// use this to temporarily copy c-blocks out of z
	Matrix<Type, _n, 1> x_bar_tmp;	// use this to temporarily copy x_bar-blocks out of z
	Matrix<Type, _nSt, 1> checkX;		// length of check vector for states
	Matrix<Type, _nInp, 1> checkU;		// vector to check input
	Matrix<Type, _nF_xTheta, 1> checkTheta;	// vector to check last constraints
	
	
	// Variables needed to represent Phi -> see documentation
	Matrix<Type, _m, _m> Omicron[_N+1];
	Matrix<Type, _n, _n> Rho[_N];
	Matrix<Type, _n, _m> Sigma[_N];
	
	// Variable to represent Y = C*Phi*C'
	Matrix<Type, _n, _n> Y[3][2*_N];	// Y[i][j] = Y_{i+1,j+1}, i=0,1,2
										// Y[1][0], Y[2][0], Y[2][1] not used
										
	// !!!!!!! Variables to represent L: L*L'=Y
	// !!!!!!! need explicitly: double, Dynamic, array-length, e.g. [2*Horizon]
	LLT<Matrix<double, Dynamic, Dynamic>, Lower> L_diag[2*50];	// ****** diagonal matrices of L computed using Cholesky
	Matrix<Type, _n, _n> L_diag_transp[2*_N];
	Matrix<Type, _n, _n> L_offDiag[2][2*_N];	// off-diag matrices are general square matrices
	Matrix<Type, _n, _n> L_offDiag_transp[2][2*_N];	// transpose of above
	
	// !!!!!!!! Variables to represent Phi = L*L'
	// !!!!!!!! need explicitly: double, Dynamic, array-length, e.g. horizon = 50
	LLT<Matrix<double, Dynamic, Dynamic>, Lower> LOmicron_diag[50+1];	//*** [horizon +1]
	Matrix<double, _m, _m> LOmicron_diag_transp[_N+1];
	LLT<Matrix<double, Dynamic, Dynamic>, Lower> LPi_diag[50];	// ***** [horizon]
	Matrix<double, _n, _n> LPi_diag_transp[_N];
	LLT<Matrix<double, Dynamic, Dynamic>, Lower> LRho_diag[50];	// ***** [horizon] 
	Matrix<double, _n, _n> LRho_diag_transp[_N];
	Matrix<Type, _m, _n> LSigma_offDiag[_N];	// last element used only if _pos_omega == _N
	Matrix<Type, _n, _m> LSigma_offDiag_transp[_N];	// last element used only if _pos_omega == _N
	Matrix<Type, _m, _n> LLambda0;		// used for LLT decomposition if _pos_omega != _N
	Matrix<Type, _n, _m> LLambda0_transp;
	Matrix<Type, _m, _m> LLambda1;		// used for LLT decomposition if _pos_omega != _N
	Matrix<Type, _m, _m> LLambda1_transp;
	
	// Variables needed to compute Y = C'*Phi_tilde * C;
	Matrix<Type, _m, _n> U[3+(_N-1)*3];		// L*U = C'; last element only used in case _pos_omega == _N
	Matrix<Type, _n, _n> U_bar[2+(_N-1)*5];	// L*U = C'
	Matrix<Type, _m, _n> UO[3];		// off-diag elements for _pos_omega != _N
	Matrix<Type, _m, _n> X[3+(_N-1)*3];		// L'*X = U
	Matrix<Type, _n, _n> X_bar[2+(_N-1)*5];	// L'*X = U
	Matrix<Type, _m, _n> XO[3];		// off-diag elements for _pos_omega != _N
	
	// residua: r_d = 2*H*z + g + kappa*P'*d + C'*nu
	//			r_p = C*z - b
	Matrix<Type, 2*_N*_n, 1> r_p;	// primal residual variable
	//Matrix<Type, 2*_N*_n, 1> b;	// RHS, needed to compute primal residuum
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> r_d;	// dual residual variable
	Matrix<Type, _N*(_nSt + _nInp) + _nF_xTheta, 1> d;	// (d)_i = 1/(h-P'z)_i
	
	// search directions dz and dnu and variables z and nu
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> z;	// z(+) = z + s*dz
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> dz;	
	Matrix<Type, 2*_N*_n, 1> nu;	//nu(+) = nu + s*dnu
	Matrix<Type, 2*_N*_n, 1> dnu;
	
	Matrix<Type, 2*_N*_n, 1> beta;	// beta = -r_p + C*Phi_tilde*r_d
	
	// pointers might be dangerous b/c no control over length
 	const Matrix<Type, _n, 1> *x_star;		// at each time step, we would like to track x_star[]
	Matrix<Type, _m, 1> u_star[_N];		// u_star[] is the least squares problem to track x_star[], but might not be feasible
	
	double reg;		// regularization term
	
	
// ----------- optimal result --------------
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> z_opt;	// contains optimal variable
	double compTime;	// stores computational time
	
// ---------- vvv  Matrices for PhaseI() -------------
	double kappa_start_PhaseI;

	// Matrix<Type, Dynamic, 1, 0, _N*(_nInp+_nSt)+_nF_xTheta, 1> gamma;
	Matrix<Type, _N*(_nInp+_nSt)+_nF_xTheta, 1> gamma;
	// Matrix<Type, Dynamic, 1, 0, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> z_hat;	// long vector
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> z_hat;	// long vector
	// Matrix<Type, Dynamic, 1, 0, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> dz_hat;	// long vector
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> dz_hat;	// long vector
	// Matrix<Type, Dynamic, 1, 0, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> r_d_hat;	// long vector
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> r_d_hat;	// long vector
	double reg_hat;		// regularisation term in Matrix H_hat
	
	// Matrix<Type, Dynamic, Dynamic, 0, _nSt, _n+_nSt> Fx_hat[_N];	// array of full-rank state constraint matrices
	Matrix<Type, _nSt, _n+_nSt> Fx_hat[_N];	// array of full-rank state constraint matrices
	// Matrix<Type, Dynamic, Dynamic, 0, _n+_nSt, _nSt> Fx_hat_transp[_N];	// transpose of above
	Matrix<Type, _n+_nSt, _nSt> Fx_hat_transp[_N];	// transpose of above
	// Matrix<Type, Dynamic, Dynamic, 0, _nInp, _m+_nInp> Fu_hat[_N];	// array of (full-rank) input constraint matrices
	Matrix<Type, _nInp, _m+_nInp> Fu_hat[_N];	// array of (full-rank) input constraint matrices
	// Matrix<Type, Dynamic, Dynamic, 0, _m+_nInp, _nInp> Fu_hat_transp[_N];	//transpose of above
	Matrix<Type, _m+_nInp, _nInp> Fu_hat_transp[_N];	//transpose of above
	Matrix<Type, _nF_xTheta, _n+_nSt> F_xTheta_hat;	// F_xTheta*x[m+N]+F_theta*theta <= f_xTheta
	Matrix<Type, _n+_nSt, _nF_xTheta> F_xTheta_hat_transp;	// transpose of above
	Matrix<Type, _nF_xTheta, _m+_nF_xTheta> F_theta_hat;	// full-rank constraint matrix on theta
	Matrix<Type, _m+_nF_xTheta, _nF_xTheta> F_theta_hat_transp;	// transpose of above
	// Matrix<Type, Dynamic, Dynamic, 0, _nInp, _n+_nSt> Fu_bar_hat[_N-1];	// array of (full-rank) input constraint matrices, #row may change, but <= nInp
	Matrix<Type, _nInp, _n+_nSt> Fu_bar_hat[_N-1];	// array of (full-rank) input constraint matrices, #row may change, but <= nInp
	// Matrix<Type, Dynamic, Dynamic, 0, _n+_nSt, _nInp> Fu_bar_hat_transp[_N-1];	//transpose of above
	Matrix<Type, _n+_nSt, _nInp> Fu_bar_hat_transp[_N-1];	//transpose of above
	
	// Matrix<Type, Dynamic, 1, 0, _m+_nInp, 1> g_hat_c[_N];	// elements in cost g_hat 
	Matrix<Type, _m+_nInp, 1> g_hat_c[_N];	// elements in cost g_hat 
	// Matrix<Type, Dynamic, 1, 0, _n+_nSt , 1> g_hat_x[_N];	// elements in cost g_hat
	Matrix<Type, _n+_nSt , 1> g_hat_x[_N];	// elements in cost g_hat
	Matrix<Type, _m+_nF_xTheta, 1> g_hat_theta;
	
	// Matrix<Type, _n, Dynamic, 0, _n, _m+_nInp> Bm_hat[_N];	// Dynamic, b/c #columns depend on # constraints
	Matrix<Type, _n, _m+_nInp> Bm_hat[_N];	// Dynamic, b/c #columns depend on # constraints
	// Matrix<Type, Dynamic, _n, 0, _m+_nInp, _n> Bm_hat_transp[_N];	// Dynamic, b/c #columns depend on # constraints
	Matrix<Type, _m+_nInp, _n> Bm_hat_transp[_N];	// Dynamic, b/c #columns depend on # constraints
	// Matrix<Type, _n, Dynamic, 0, _n, _m+_nInp> B_hat[_N];
	Matrix<Type, _n, _m+_nInp> B_hat[_N];
	// Matrix<Type, Dynamic, _n, 0, _m+_nInp, _n> B_hat_transp[_N];
	Matrix<Type, _m+_nInp, _n> B_hat_transp[_N];
	// Matrix<Type, _n, Dynamic, 0, _n, _n+_nSt> Bm_bar_hat[_N-1];
	Matrix<Type, _n, _n+_nSt> Bm_bar_hat[_N-1];
	// Matrix<Type, Dynamic, _n, 0, _n+_nSt, _n> Bm_bar_hat_transp[_N-1];
	Matrix<Type, _n+_nSt, _n> Bm_bar_hat_transp[_N-1];
	// Matrix<Type, _n, Dynamic, 0, _n, _n+_nSt> A_bar_hat[_N-1];
	Matrix<Type, _n, _n+_nSt> A_bar_hat[_N-1];
	// Matrix<Type, Dynamic, _n, 0, _n+_nSt, _n> A_bar_hat_transp[_N-1];
	Matrix<Type, _n+_nSt, _n> A_bar_hat_transp[_N-1];
	// Matrix<Type, _n, Dynamic, 0, _n, _n+_nSt> Identity_hat[_N];
	Matrix<Type, _n, _n+_nSt> Identity_hat[_N];
	// Matrix<Type, Dynamic, _n, 0, _n+_nSt, _n> Identity_hat_transp[_N];
	Matrix<Type, _n+_nSt, _n> Identity_hat_transp[_N];
	
	Matrix<Type, Dynamic, Dynamic> Omicron_hat[_N+1];	// upper bound either _m+_nInp or _m+_nF_xTheta
	// Matrix<Type, Dynamic, Dynamic, 0, _n+_nSt, _n+_nSt> Rho_hat[_N];
	Matrix<Type, _n+_nSt, _n+_nSt> Rho_hat[_N];
	Matrix<Type, Dynamic, Dynamic> Sigma_hat[_N];		// upper bound either ...+_nInp or ...+_nF_xTheta
	
	Matrix<Type, Dynamic, Dynamic> LOmicron_hat_diag_transp[_N+1];
	Matrix<Type, _n, _n> LPi_hat_diag_transp[_N];	// stores the transpose
	Matrix<Type, _n+_nSt, _n+_nSt> LRho_hat_diag_transp[_N];
	// Matrix<Type, Dynamic, Dynamic, 0, _n+_nSt, _n+_nSt> LRho_hat_diag_transp[_N];
	Matrix<Type, Dynamic, Dynamic> LSigma_hat_offDiag[_N];
	Matrix<Type, Dynamic, Dynamic> LSigma_hat_offDiag_transp[_N];
	Matrix<Type, Dynamic, Dynamic> LLambda0_hat;		// used for LLT decomposition if _pos_omega != _N
	Matrix<Type, Dynamic, Dynamic> LLambda0_hat_transp;
	Matrix<Type, Dynamic, Dynamic> LLambda1_hat;		// used for LLT decomposition if _pos_omega != _N
	Matrix<Type, Dynamic, Dynamic> LLambda1_hat_transp;
	
	Matrix<Type, Dynamic, _n> U_hat[3+(_N-1)*3];		// L*U = C', needed to compute Y
	Matrix<Type, Dynamic, _n> U_bar_hat[2+(_N-1)*5];	// L*U = C'
	Matrix<Type, Dynamic, _n> UO_hat[3];
	Matrix<Type, Dynamic, _n> X_hat[3+(_N-1)*3];		// L'*X = U
	Matrix<Type, Dynamic, _n> X_bar_hat[2+(_N-1)*5];	// L'*X = U
	Matrix<Type, Dynamic, _n> XO_hat[3];
	
	
	double difference;
// ------------ ^^^ end Matrices PhaseI() --------------	
	
	
	
	
	
	
	

// ---------- private methods ----------
	// "testPos" tests if the given z_warm satisfies P*z_warm < h								
	bool testPos();		// returns 1, inequality is satisfied
						// return 0, if inequality is violated
	
	void phaseI();	// in case P*z_warm >= h, phaseI will retrieve a good starting point
	void compD();	// compute "d"-vector
	void compRdRp(); 	// compute the dual and primal residua
	void compDzDnu();
	void compPhi();
	void compPhi_tilde();
	void compY();
	void compBeta();		// beta = -r_p + C*Phi_tilde*r_d;
	void compL();		// L*L' = Y
	void compDnu();		// Y*Dnu = 
	void compDz();
	double compT();
	double getTimer();	// method returning the computational time
	void compRQ();		// computes the cost vectors r and q

// ------- vvv methods for PhaseI() ----------
	void compMatrices_PhaseI();
	void compGamma_PhaseI();
	void compRdRp_PhaseI();
	void compDz_hatDnu_hat_PhaseI();
	void compPhi_hat_PhaseI();
	void compPhi_hat_tilde_PhaseI();
	void compY_hat_PhaseI();
	void compBeta_hat_PhaseI();
	void compDz_hat_PhaseI();
	bool testPos_PhaseI();
	void compD_hat_PhaseI();
	bool testNeg_PhaseI();
	double compT_PhaseI();		// computes largest t: P_hat*(z_hat+t*dz_hat) < h
// --------- ^^^ functions for PhaseI() -------


	
public:
	LBmpcTP();	// default constructor
	
	// usual constructor: uses references whenever possible to avoid copying of big matrices
	// note: arrays cannot be copied by value, always given as pointers, length implicitly given by _N
	LBmpcTP(double kappa_arg, double kappa_PhaseI_arg, int n_iter_arg, double mu_arg, double eps_barrier_arg, double eps_nt_arg, double eps_normRp_arg,
		  	double eps_ls_arg, double alpha_ls_arg, double beta_ls_arg, double reg_arg, double reg_PhaseI_arg,
		   const Matrix<Type, _n, _n> &A_arg, Matrix<Type, _n, _m> &B_arg,
		   const Matrix<Type, _n, _n> &Q_tilde_arg, const Matrix<Type, _n, _n> &Q_tilde_f_arg,
		   const Matrix<Type, _m, _m> &R_arg, const Matrix<Type, _nSt, _n> Fx_arg[], 
		   const Matrix<Type, _nSt, 1> fx_arg[], const Matrix<Type, _nInp, _m> Fu_arg[],
		   const Matrix<Type, _nInp, 1> fu_arg[], const Matrix<Type, _nF_xTheta, _n> &F_xTheta_arg,
		   const Matrix<Type, _nF_xTheta, _m> &F_theta_arg, const Matrix<Type, _nF_xTheta, 1> &f_xTheta_arg,
		   const Matrix<Type, _m, _n> &K_arg, const Matrix<Type, _n, 1> &s_arg
	      );	
	

	// "step" computes and returns the optimal z-vector containing the optimal offset c[m]
	Matrix<Type, _m, 1> step(const Matrix<Type, _n, _n> &Lm_arg, const Matrix<Type, _n, _m> &Mm_arg, const Matrix<Type, _n, 1> &tm_arg,
														const Matrix<Type, _n, 1> &x_hat_arg, Matrix<Type, _N*(_m + _n + _n) + _m, 1> &z_warm_arg,
														const Matrix<Type, _n, 1> x_star_arg[]
													   );
													
	//	const Matrix<Type, _n, 1> &q_tilde_arg, const Matrix<Type, _n, 1> &q_tilde_f_arg,
	// const Matrix<Type, _m, 1> &r_arg,
	
	//~LBmpcTP();	// destructor
};




//  ==================== Implementation of Methods ==================

template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::LBmpcTP()
{
	// do nothing
}

template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::LBmpcTP(double kappa_arg, double kappa_PhaseI_arg, int n_iter_arg, 
	   double mu_arg, double eps_barrier_arg, double eps_nt_arg, double eps_normRp_arg, double eps_ls_arg, double alpha_ls_arg, double beta_ls_arg,
	   double reg_arg, double reg_PhaseI_arg, const Matrix<Type, _n, _n> &A_arg, Matrix<Type, _n, _m> &B_arg,
	   const Matrix<Type, _n, _n> &Q_tilde_arg, const Matrix<Type, _n, _n> &Q_tilde_f_arg,
	   const Matrix<Type, _m, _m> &R_arg, const Matrix<Type, _nSt, _n> Fx_arg[], 
	   const Matrix<Type, _nSt, 1> fx_arg[], const Matrix<Type, _nInp, _m> Fu_arg[],
	   const Matrix<Type, _nInp, 1> fu_arg[], const Matrix<Type, _nF_xTheta, _n> &F_xTheta_arg,
	   const Matrix<Type, _nF_xTheta, _m> &F_theta_arg, const Matrix<Type, _nF_xTheta, 1> &f_xTheta_arg,
	   const Matrix<Type, _m, _n> &K_arg, const Matrix<Type, _n, 1> &s_arg
      )
{
	
	kappa_start = kappa_arg;
	kappa_start_PhaseI = kappa_PhaseI_arg;
	n_iter = n_iter_arg;
	mu = mu_arg;
	eps_barrier = eps_barrier_arg;
	eps_nt = eps_nt_arg;
	eps_ntSq = eps_nt*eps_nt;
	eps_normRp = eps_normRp_arg;
	eps_normRpSq = eps_normRp*eps_normRp;
	eps_ls = eps_ls_arg;
	reg = reg_arg;
	reg_hat = reg_PhaseI_arg;
	t = 1;
	alpha_ls = alpha_ls_arg;
	beta_ls = beta_ls_arg;
	offset = _n + _n + _m;
	
	
	
	A = A_arg;
	A_transp = A_arg.transpose();
	B = B_arg;
	B_transp = B_arg.transpose();
	Q_tilde = Q_tilde_arg;
	Q_tilde_f = Q_tilde_f_arg;
	R = R_arg;	
	K = K_arg;
	K_transp = K.transpose();
		
	// count number of constraints and save array
	// num_constr = 0;
	for (int i=0; i < _N; i++)
	{
		Fx[i] = Fx_arg[i];
		Fx_transp[i] = Fx_arg[i].transpose();
		fx[i] = fx_arg[i];
		// Fx_rows[i] = fx[i].rows();
		
		// num_constr = num_constr + Fx_rows[i];
		
		Fu[i] = Fu_arg[i];
		Fu_transp[i] = Fu_arg[i].transpose();
		fu[i] = fu_arg[i];
		// Fu_rows[i] = fu[i].rows();
		Fu_bar[i] = Fu[i]*K;		// the first one is NOT used
		Fu_bar_transp[i] = Fu_bar[i].transpose();
		
		// num_constr = num_constr + Fu_rows[i];
	}
	
	F_xTheta = F_xTheta_arg;
	F_xTheta_transp = F_xTheta_arg.transpose();
	F_theta = F_theta_arg;
	F_theta_transp = F_theta_arg.transpose();
	f_xTheta = f_xTheta_arg;
	// num_constr = num_constr + _nF_xTheta;
	
	num_constr = _N*(_nInp + _nSt) + _nF_xTheta;
	
	s = s_arg;

	// d.resize(_N*(_nInp + _nSt) + _nF_xTheta);
	
	//-------- do some computations	
	S = K.transpose() * R;
	S_transp = S.transpose();
	Q_bar = S * K;
	// q_bar = K.transpose() * r;
	A_bar = A + B*K;
	A_bar_transp = A_bar.transpose();	
}

// step function returns optimal input
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
Matrix<Type, _m, 1> LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::step(const Matrix<Type, _n, _n> &Lm_arg, const Matrix<Type, _n, _m> &Mm_arg, const Matrix<Type, _n, 1> &tm_arg,
																									   const Matrix<Type, _n, 1> &x_hat_arg, Matrix<Type, _N*(_m + _n + _n) + _m, 1> &z_warm_arg,
																									  const Matrix<Type, _n, 1> x_star_arg[]
												   													  )
// const Matrix<Type, _n, 1> &q_tilde_arg, const Matrix<Type, _n, 1> &q_tilde_f_arg, const Matrix<Type, _m, 1> &r_arg
{
	// initialization
	Lm = Lm_arg;
	Mm = Mm_arg;
	tm = tm_arg;
	x_hat = x_hat_arg;
	z_warm = z_warm_arg;
	kappa = kappa_start;
	x_star = x_star_arg;
	tm_tilde = tm_arg + s;
	Am_tilde = A + Lm;
	Am_tilde_transp = Am_tilde.transpose();
	Bm = B + Mm;
	Bm_transp = Bm.transpose();
	Bm_bar = Bm * K;
	Bm_bar_transp = Bm_bar.transpose();
	
	compRQ();	// compute u_star, x_star -> cost matrices 
	// test, if given w_warm satisfies (P*z_warm<h), otherwise better z_warm is computed using phaseI()
	// testPos(): returns 1 if inequality is satisfied
	//					  0 if inequality not satisfied
	z = z_warm;
	// cout << "z_warm:" << endl << z_warm << endl << endl;
	
	nu.setConstant(1);	// any nu good to initialize
	srand((unsigned)time(0));
	nu.setRandom();
	nu = 100*nu;
	
	if( !testPos() )	// testPos() operates on z
	{
			phaseI();		// sets new z_warm and an appropriate nu
			// return K*x_hat + z.template segment<_m>(0);
	}
	else
	{
		cout << "chosen z_warm is in domain." << endl;
	}
	// return K*x_hat + z.template segment<_m>(0);

	kappa = kappa_start;
	z = z_warm;		// the (new) z_warm satisfies inequality
	compD();	// compute d-vector
	// some auxiliary variables
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> z_orig;	// stores the "old z", z stores z = z_tmp + t*dz
	Matrix<Type, 2*_N*_n, 1> nu_orig;
	double resNormSq_orig;
	double resNormSq;	// squared
	double resNormSqRp;
	bool isFeasible;	// stores if P*(z+t*dz) < h
	bool cond;
	int i;		// counter for backtracking line search
	int j;		// number of newton iterations for a fixed kappa
	int loopCounter = 0;
	int i_kappa = 0;	// delete it later
	do 			// kappa changes
	{
		i_kappa++;	
		cout << "============= kappa: " << kappa << ", i_kappa: " << i_kappa << " ================" << endl;
		compRdRp();		// z is already in domain, using new kappa to compute primal and dual residua
		j = 0;
		do 		// kappa fixed
		{
			j++;	// number of Newton steps for fixed kappa
			// cout << "----------- kappa fixed, round: " << j << "--------------" << endl << endl;			
			compDzDnu();
			
			// --------- line search part --------------
			// get smaller and smaller t until following two conditions are satisfied
			// 1) ||r_d(z+t*dz,nu+t*dnu)||_2 + ||r_p(z+t*dz,nu+t*dnu)||_2 <= (1-alpha*t)*(||r_d||_2 + ||r_p||_2)
			// 2) P*(z+t*dz) < h
			t = compT();	// retrieve largest t s.t. P*(z+t*dz) < h			
			if (t < 0)
			{
				cout << "! compT() returned negative t:" << t << endl << endl;
				// cout << "z_hat:" << endl << z_hat << endl << endl;
				// cout << "dz_hat:" << endl << dz_hat << endl << endl;
				return  K*x_hat + z.template segment<_m>(0);
			}
			// t = 1;
			z_orig = z;	// stores the "old z", z stores z = z_tmp + t*dz
			nu_orig = nu;
			
			resNormSq_orig = r_d.squaredNorm() + r_p.squaredNorm();
			cond = 1;
			i = 0;
			loopCounter++;
			isFeasible = 0;
			
			while (cond)
			{
				// loopCounter++;
				i++;
				// cout << "entering round: " << i << endl << endl;
				z = z_orig + t*dz;
				// cout << "updated z with t = " << t << ":" << endl << z << endl << endl;
				isFeasible = isFeasible || testPos();		// test positivity using the new z
				if(!isFeasible)
				{
					cout << "!!! not feasible despite compT()!" << endl;
				}
				// cout << "cond1: " << cond << endl << endl;
				if (isFeasible)	// means that z is in domain, only then check if inequality holds
				{
					nu = nu_orig + t*dnu;
					compD();
					compRdRp();		// class variables r_d and r_p are updated
					resNormSqRp = r_p.squaredNorm();
					resNormSq = r_d.squaredNorm() + resNormSqRp;
					cond = ( (resNormSq > (1-alpha_ls*t)*(1-alpha_ls*t)*resNormSq_orig) );
				}		
				t = beta_ls*t;
				if (t <= eps_ls)
				{
					break;	// stop, no improvement
				}
			}	// z and nu have automatically been updated
		
			cout << "line search needed " << i << " steps with step size t: " << t/beta_ls << endl << endl;
			
			if (j >= n_iter && i_kappa > 0)		// i_kappa==1 might need more steps
			{
				cout << "!!!!!! more than " << n_iter <<" steps required !!!!!!" << endl;
				cout << "loopCounter: " << loopCounter << endl;
				break;
			}
		} while( (resNormSq > eps_ntSq) || (resNormSqRp > eps_normRpSq) );
		// while( (resNormSq > eps_ntSq) || (resNormSqRp > eps_normRpSq) );
		cout << "needed " << j << " rounds" << endl << endl;
		kappa = kappa*mu;
	} while( kappa/mu*num_constr > eps_barrier );
	//while( kappa/mu*num_constr > eps_barrier );
	
	cout << " =====> computed optimal z_vector:" << endl << setprecision (15) << z << endl << endl;
	cout << setprecision (15) << "kappa: " << kappa << " after i_kappa: " << i_kappa << endl;
	cout << setprecision (15) << "eps_barrier: " << eps_barrier << endl;
	cout << setprecision (15) << "kappa/mu*num_constr: " << kappa/mu*num_constr << endl << endl;
	cout << setprecision (15) << "loopCounter: " << loopCounter << endl << endl;
	
	return K*x_hat + z.template segment<_m>(0);
}


// ------------ function tests if z_warm satisfies P*z < h --------
// ------- might consider computing entire matrix check and then iterate for neg. elements ----
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: testPos()
{
	
	// cout << endl << "======= starting testPos() =======" << endl;
	
	// special treatment at beginning
	c_tmp = z.template segment<_m>(0);
	checkU = (fu[0] - Fu_bar[0]*x_hat) - (Fu[0]*c_tmp);	// should be >0
	
	for (int j=1; j <= _nInp; j++)
	{
		if (checkU[j-1] <= 0)
		{
			// cout << "fu Problem at beginning at" << j << endl << endl;
			return 0;
		}
	}
	// cout << "finished special treatment at start" << endl;
	
	// class variable: offset = _n + _n + _m;
	// general treatment in the middle 
	for (int i=1; i <= _N-1; i++) 	// blocks in the middle
	{	// compute (h - P*z)_i
		// cout << "round " << i << "in general treatment" << endl;
		x_bar_tmp = z.template segment<_n>((i-1)*offset+_m+_n);
		c_tmp = z.template segment<_m>((i-1)*offset+_m+_n+_n);
		checkX = fx[i-1] - Fx[i-1]*x_bar_tmp;
		// cout << "computed check" << endl << check << endl << endl;
		for (int j=1; j<= _nSt; j++)
		{
			if (checkX[j-1] <= 0)
			{
				// cout << "fx Problem at" << i << "at " << j << endl << endl;
				return 0;
			}
		}
		
		checkU = fu[i] - (Fu_bar[i]*x_bar_tmp + Fu[i]*c_tmp);
		for (int j=1; j<=_nInp; j++)
		{
			if (checkU[j-1] <= 0)
			{
				return 0;
			}
		}
		//cout << "for-loop completion of round" << i << endl;
	}
	
	// cout << "finished general treatment" << endl;
	
	// special case for last blocks
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);	// depends on where it is
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	checkX = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	// cout << check << endl << endl;
	for (int j=1; j<= _nSt; j++)
	{
		if (checkX[j-1] <= 0)
		{
			return 0;
		}
	}
	
	// cout << "leihou" << endl;
	
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);	// depends on position of invariant set
	checkTheta = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	// cout << check << endl << endl;
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if (checkTheta[j-1] <= 0)
		{
			// cout << "f_xTheta Problem in last block at " << j << endl << endl;
			return 0;
		}
	}
	// cout << "======= leaving testPos() ========" << endl << endl;
	
	return 1;	// no problems detected
}


// ------------ function finds z_warm such that P*z_warm < h -----------
// ------------ implementing basically the other algorithms ------------
// solve the problem: min{1'*gamma + z'*H_hat*z}, s.t. {P*z-h < gamma, C*z = b}
// reformulate as: min{z_hat*H_hat*z + g_hat'*z_hat}, s.t. {P_hat*z_hat < h, C_hat * z_hat = b}
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: phaseI()
{
	cout << "---------- Pz < h not satisfied, starting phaseI() -----------" << endl;
	/* definition of all the variables needed */
	// z_hat.resize(num_constr + _N*(_m + _n + _n) + _m);
	// dz_hat.resize(num_constr + _N*(_m + _n + _n) + _m);
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> z_hat_orig;	// vector storing the original data
	// z_hat_orig.resize(num_constr + _N*(_m + _n + _n) + _m);
	// r_d_hat.resize(num_constr + _N*(_m + _n + _n) + _m);
	Matrix<Type, 2*_N*_n, 1> nu_orig;
	// gamma.resize(num_constr);
	
	nu.setConstant(1);	// any nu good to initialize
	// nu << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
	nu << 8, 1, 7, 8, 5, 0, 4, 10, 1, 8, 4, 3, 5, 4, 9, 4, 1, 3, 2, 8, 6, 7,
	   		7, 5, 3, 1, 8, 3, 6, 5, 7, 3, 3, 7, 10, 3, 9, 9, 8, 8, 3, 9, 8, 8,
			8, 5, 9, 9, 0, 7, 5, 6, 10, 8, 10, 3, 3, 3, 3, 5, 3, 9, 7, 8, 2, 1,
			1, 6, 6, 4, 0, 2, 5, 8, 8, 9, 8, 5, 3, 1, 6, 2, 1, 8, 5, 7, 1, 4,
			6, 9, 9, 4, 4, 1, 7, 2, 9, 4, 5, 2;
	// srand((unsigned)time(0));
	// nu.setRandom();
	// nu = 10*nu;
	
	compMatrices_PhaseI();	// computes Fx_hat, Fu_hat, Fu_bar_hat, g_hat_c, g_hat_x, Bm_hat, Identity_hat...

	difference = 100;	// h-P*z_warm = difference
	compGamma_PhaseI();
	// cout << setprecision(15) << "gamma:" << endl << gamma << endl << endl;
	// return;
	
	/* implementation of functions */
	d.setConstant(1/difference);	// 1/difference doesn't work, returns 0....; difference set in compGamma_PhaseI
	
	// vvv build z_hat
	int tmp = 0;	// counter for position 
	for (int i=1; i<= _N; i++)
	{
		z_hat.template segment<_m>( (i-1)*offset+tmp ) = z.template segment<_m>( (i-1)*offset );
		z_hat.template segment<_nInp>((i-1)*offset+_m+tmp) = gamma.template segment<_nInp>(tmp);
		tmp = tmp + _nInp;
		z_hat.template segment<_n>( (i-1)*offset+tmp+_m ) = z.template segment<_n>( (i-1)*offset+_m );
		z_hat.template segment<_n>( (i-1)*offset+tmp+_m+_n ) = z.template segment<_n>((i-1)*offset+_m+_n);
		z_hat.template segment<_nSt>((i-1)*offset+tmp+_m+_n+_n) = gamma.template segment<_nSt>(tmp);
		tmp = tmp + _nSt;
	}
	z_hat.template segment<_m>(_N*offset+tmp) = z.template segment<_m>(_N*offset);
	z_hat.template segment<_nF_xTheta>(_N*offset+_m+tmp) = gamma.template segment<_nF_xTheta>(tmp);
	// cout << setprecision(15) << "z_hat" << endl << z_hat << endl << endl;
	// return;
	// ^^^ computed z_hat
		
	double resNormSq_orig;
	double resNormSq;
	double resNormSqRp;
	int j;
	int loopCounter = 0;
	int i_kappa = 0;
	kappa = kappa_start_PhaseI;
	bool cond;
	bool isFeasible;
	int i;
	do
	{
		i_kappa++;
		cout << "============= kappa PhaseI: " << kappa << ", i_kappa: " << i_kappa << " ================" << endl;
		compRdRp_PhaseI();		// computes r_p and r_d_hat
		// cout << setprecision(15) << "r_p in PhaseI()" << endl << r_p << endl << endl;
		// cout << "r_d_hat in PhaseI()" << endl << r_d_hat << endl << endl;
		// return;
		
		j = 0;
		do  	// kappa fixed 
		{
			j++;
			// cout << "----------- kappa fixed PhaseI(), round: " << j << "--------------" << endl << endl;			
			compDz_hatDnu_hat_PhaseI();	// compute dz_hat and dnu
			// cout << setprecision(15) << "dz_hat" << endl << dz_hat << endl << endl;
			// cout << setprecision(15) << "dnu" << endl << dnu << endl << endl;
			// return;
			
			// ---- line search part --------
			// decrease t until both conditions are satisfied
			// 1) ||res(z+t*dz,nu+t*dnu)||_2 <= (1-alpha*t)*(||r_d||_2 + ||r_p||_2)
			// 2) P*(z+t*dz) < h
			// t = 1;
			// cout << "before compT_PhaseI" << endl << endl;
			t = compT_PhaseI();		// computes largest t: P_hat*(z_hat+t*dz_hat) < h
			// cout << setprecision(15) << "z_hat" << endl << z_hat << endl << endl;
			// cout << "dz_hat" << endl << dz_hat << endl << endl;
			// cout << setprecision(15) << "compT_PhaseI: " << t << endl << endl;
			// return;
			if (t < 0)
			{
				cout << "compT_PhaseI() returned negative t:" << t << endl << endl;
				cout << "z_hat:" << endl << z_hat << endl << endl;
				cout << "dz_hat:" << endl << dz_hat << endl << endl;
				return;
			}
			z_hat_orig = z_hat;
			nu_orig = nu;
			resNormSq_orig = r_d_hat.squaredNorm() + r_p.squaredNorm();

			cond = 1;
			isFeasible = 0;
			i = 0;
			loopCounter++;
			while (cond)
			{
				// loopCounter++;
				i++;
				z_hat = z_hat_orig + t*dz_hat;	// guaranteed to be feasible, b/c t is small enough
				isFeasible = isFeasible || testPos_PhaseI();
				if(!isFeasible)
				{
					cout << "not feasible!" << endl;
					cout << "z_hat_orig:" << endl << z_hat_orig << endl << endl;
					cout << "t: " << t << endl;
					cout << "dz_hat" << endl << dz_hat << endl << endl;
					cout << "z_hat:" << endl << z_hat << endl << endl;
					return;
				}
				// cout << "cond at round " << i << " is: " << cond << endl << endl;
				if (isFeasible)
				{
					// stop if gamma < 0 and z_warm is computed automatically
					if( testNeg_PhaseI() )
					{
						cout << setprecision(15) << "new z_warm is: " << endl << z_warm << endl << endl;
						// cout << "d:" << endl << d << endl << endl;
						cout << "loopCounter needed for finding z_warm: " << loopCounter << endl << endl;
						return;
					}
					nu = nu_orig + t*dnu;
					compD_hat_PhaseI();	// stored in d 
					// cout << setprecision(15) << "d in phaseI" << endl << d << endl << endl;
					compRdRp_PhaseI();
					// cout << setprecision(15) << "r_p in PhaseI()" << endl << r_p << endl << endl;
					// cout << "r_d_hat in PhaseI()" << endl << r_d_hat << endl << endl;
					// return;
					resNormSqRp = r_p.squaredNorm();
					resNormSq = r_d_hat.squaredNorm() + resNormSqRp;
					cond = ( (resNormSq > (1-alpha_ls*t)*(1-alpha_ls*t)*resNormSq_orig) );
				}
				t = beta_ls*t;
				// cout << "t: " << t << endl << endl;
			}	
			// cout << "line search in PhaseI needed " << i << " steps with step size t: " << t/beta_ls << endl << endl;
			if (j == 200)
			{
				cout << "!!!!!! PhaseI more than 200 steps required !!!!!!" << endl;
				cout << "loopCounter: " << loopCounter << endl;
				// cout << "z_hat:" << endl << z_hat << endl << endl;
				return;
			}
		} while( (resNormSq > eps_ntSq) || (resNormSqRp > eps_normRpSq) );
		// while( (resNormSq > eps_ntSq) || (resNormSqRp > eps_normRpSq) );
		cout << "needed " << j << " rounds" << endl << endl;
		kappa = kappa*mu;
	} while( kappa/mu*num_constr > eps_barrier );
	cout << "***** phaseI() did NOT converge." << endl;
	//while( kappa/mu*num_constr > eps_barrier );
	return;
}


// ------------ function computes primal and dual residua  -------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compRdRp()
{
	//cout << "========= starting compRdRp() =========" << endl;
	
	// ---------- compute primal dual r_p = C*z - b, done two blocks
	//special treatment at the beginning
	r_p.template segment<_n>(0) = -Bm*z.template segment<_m>(0) + z.template segment<_n>(_m) - ( Am_tilde + Bm_bar)*x_hat - tm_tilde;
	r_p.template segment<_n>(_n) = -B*z.template segment<_m>(0) + z.template segment<_n>(_m+_n) - (A_bar*x_hat + s);
	
	// deal with the rest, class variable: offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++)
	{
		r_p.template segment<_n>(2*_n*i) = -Am_tilde*z.template segment<_n>((i-1)*offset+_m) - Bm_bar*z.template segment<_n>((i-1)*offset+_m+_n) -
			Bm*z.template segment<_m>((i-1)*offset+_m+_n+_n) + z.template segment<_n>((i-1)*offset+_m+_n+_n+_m) - tm_tilde;
		r_p.template segment<_n>(2*_n*i+_n) = -A_bar*z.template segment<_n>((i-1)*offset+_m+_n) - B*z.template segment<_m>((i-1)*offset+_m+_n+_n) + 
						z.template segment<_n>((i-1)*offset+_m+_n+_n+_m+_n) - s;
		
	}
	
	// cout << "r_p" << endl << r_p << endl << endl;
	// return;

	// ------------- compute dual r_d = 2*H*z + g + kappa*P'*d + C'*nu;
	// handle first case separately
	r_d.template segment<_m>(0) = 2*R*z.template segment<_m>(0) + (r_vec[0]+2*S_transp*x_hat) + 
				kappa*Fu_transp[0]*d.template segment<_nInp>(0) - Bm_transp*nu.template segment<_n>(0) - B_transp*nu.template segment<_n>(_n);
						
	// handle the cases in the middle, without the end
	// three block to deal with in each round
	int offset1 = 2*_n;	// offset required in nu for C'*nu
	// int tmp = Fu_rows[0];	// offset required in d for P'*d
	for (int i=1; i<= _N-1; i++)
	{
		r_d.template segment<_n>(_m+(i-1)*offset) = 2*Q_tilde*z.template segment<_n>(_m+(i-1)*offset) + q_tilde_vec[i-1] +
				nu.template segment<_n>((i-1)*offset1) - Am_tilde_transp*nu.template segment<_n>((i-1)*offset1+_n+_n);
		
		if (i != _pos_omega)	
		{	
			r_d.template segment<_n>(_m+(i-1)*offset+_n) = 2*Q_bar*z.template segment<_n>(_m+(i-1)*offset+_n) + 2*S*z.template segment<_m>(_m+(i-1)*offset+_n+_n) + q_bar_vec[i] +
					// kappa*Fx_transp[i-1]*d.segment(tmp,Fx_rows[i-1]) + kappa*Fu_bar_transp[i]*d.segment(tmp+Fx_rows[i-1],Fu_rows[i]) + 
					kappa*Fx_transp[i-1]*d.template segment<_nSt>(_nInp+(i-1)*(_nInp+_nSt)) + kappa*Fu_bar_transp[i]*d.template segment<_nInp>(i*(_nInp+_nSt)) + 
					nu.template segment<_n>((i-1)*offset1+_n) - Bm_bar_transp*nu.template segment<_n>((i-1)*offset1+_n+_n) - A_bar_transp*nu.template segment<_n>((i-1)*offset1+_n+_n+_n);
		}
		else	// must add the additional term: F_xTheta'*d(.)
		{
			r_d.template segment<_n>(_m+(_pos_omega-1)*offset+_n) = 2*Q_bar*z.template segment<_n>(_m+(_pos_omega-1)*offset+_n) + 2*S*z.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) + q_bar_vec[i] +
					// kappa*Fx_transp[_pos_omega-1]*d.segment(tmp,Fx_rows[_pos_omega-1]) + kappa*Fu_bar_transp[_pos_omega]*d.segment(tmp+Fx_rows[_pos_omega-1],Fu_rows[_pos_omega]) + 
					kappa*Fx_transp[_pos_omega-1]*d.template segment<_nSt>(_nInp+(i-1)*(_nInp+_nSt)) + kappa*Fu_bar_transp[_pos_omega]*d.template segment<_nInp>(i*(_nInp+_nSt)) + 
					nu.template segment<_n>((_pos_omega-1)*offset1+_n) - Bm_bar_transp*nu.template segment<_n>((_pos_omega-1)*offset1+_n+_n) - A_bar_transp*nu.template segment<_n>((_pos_omega-1)*offset1+_n+_n+_n)
					+ F_xTheta_transp * d.template segment<_nF_xTheta>(num_constr-_nF_xTheta);
		}
		
		// tmp = tmp + Fx_rows[i-1];
		r_d.template segment<_m>(_m+(i-1)*offset+_n+_n) = 2*S_transp*z.template segment<_n>(_m+(i-1)*offset+_n) + 2*R*z.template segment<_m>(_m+(i-1)*offset+_n+_n) + r_vec[i] +
				// kappa*Fu_transp[i]*d.segment(tmp,Fu_rows[i]) - 
				kappa*Fu_transp[i]*d.template segment<_nInp>(i*(_nInp+_nSt)) - 
				Bm_transp*nu.template segment<_n>((i-1)*offset1+_n+_n) - B_transp*nu.template segment<_n>((i-1)*offset1+_n+_n+_n);		
		// tmp = tmp + Fu_rows[i];
	}
	r_d.template segment<_n>(_m+(_N-1)*offset) = 2*Q_tilde_f*z.template segment<_n>(_m+(_N-1)*offset) + q_tilde_vec[_N-1] + nu.template segment<_n>((_N-1)*offset1);
	
	if(_pos_omega == _N)
	{
		// r_d.template segment<_n>(_m+(_N-1)*offset+_n) = kappa*Fx_transp[_N-1]*d.segment(tmp,Fx_rows[_N-1]) + kappa*F_xTheta_transp*d.segment(tmp+Fx_rows[_N-1],_nF_xTheta) + 
		r_d.template segment<_n>(_m+(_N-1)*offset+_n) = kappa*Fx_transp[_N-1]*d.template segment<_nSt>(_nInp+(_N-1)*(_nInp+_nSt)) + kappa*F_xTheta_transp*d.template segment<_nF_xTheta>(_N*(_nInp+_nSt)) + 
				nu.template segment<_n>((_N-1)*offset1+_n);
	}
	else	//standard
	{
		// r_d.template segment<_n>(_m+(_N-1)*offset+_n) = kappa*Fx_transp[_N-1]*d.segment(tmp,Fx_rows[_N-1]) + nu.template segment<_n>((_N-1)*offset1+_n);
		r_d.template segment<_n>(_m+(_N-1)*offset+_n) = kappa*Fx_transp[_N-1]*d.template segment<_nSt>(_nInp+(_N-1)*(_nInp+_nSt)) + nu.template segment<_n>((_N-1)*offset1+_n);
	}
	// tmp = tmp + Fx_rows[_N-1];
	// r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n) = kappa * F_theta_transp*d.template segment<_nF_xTheta>(tmp);
	r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n) = kappa * F_theta_transp*d.template segment<_nF_xTheta>(_N*(_nSt+_nInp));
	
	// cout << "r_p:" << endl << r_p << endl << endl;
	// cout << "r_d" << endl << r_d << endl << endl;
}


// ------------ function computes vector (d)_i=(h-P*z)_i  --------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compD()
{
	// Matrix<Type, _m, 1> c_tmp;		// use this to temporarily copy c-blocks out of z
	// Matrix<Type, _n, 1> x_bar_tmp;	// use this to temporarily copy x_bar-blocks out of z
	// Matrix<Type, Dynamic, 1> diff;		// length of vector not yet fixed
	// int tmp = 0;		// variable used to count position

	// special treatment at beginning
	c_tmp = z.template segment<_m>(0);
	checkU = (fu[0] - Fu[0]*K*x_hat) - (Fu[0]*c_tmp);	// should be >0
	for (int j=1; j <= _nInp; j++)
	{
		d[j-1] = 1/checkU[j-1];
	}
	// tmp = tmp + Fu_rows[0];
	// general treatment in the middle, class variable: offset = _n + _n + _m
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	// compute (h - P*z)_i
		x_bar_tmp = z.template segment<_n>((i-1)*offset+_m+_n);
		c_tmp = z.template segment<_m>((i-1)*offset+_m+_n+_n);		
		checkX = fx[i-1] - Fx[i-1]*x_bar_tmp;
		for (int j=1; j<= _nSt; j++)
		{
			d[_nInp+(i-1)*(_nInp+_nSt)+j-1] = 1/checkX[j-1];
		}
		// tmp = tmp + Fx_rows[i-1];
		
		checkU = fu[i] - (Fu[i]*K*x_bar_tmp + Fu[i]*c_tmp);
		for (int j=1; j<=_nInp; j++)
		{
			d[_nInp+(i-1)*(_nInp+_nSt)+_nSt+j-1] = 1/checkU[j-1];
		}
		// tmp = tmp + Fu_rows[i];
	}
	
	// special case for last blocks
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	
	checkX = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	for (int j=1; j<= _nSt; j++)
	{
		d[_nInp+(_N-1)*(_nSt+_nInp)+j-1] = 1/checkX[j-1];
	}
	// tmp = tmp + Fx_rows[_N-1];
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);
	checkTheta = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		d[_N*(_nSt+_nInp)+j-1] = 1/checkTheta[j-1];
	}
	// cout << "d: " << endl << d << endl << endl;
}

// ------- function computes primal and dual Newton steps -------------
// ------------ implementing ------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDzDnu()
{
	compPhi();	// Phi = 2*H + kappa*P'*diag(d)^2*P
	compPhi_tilde();	// Phi_tilde = Phi^{-1}; computes Chol-Dec
	compY();			// Y = C*Phi_tilde*C'
	compBeta();		// beta = -r_p + C*Phi_tilde*r_d;
	//    ---------------- compute elements of L-matrix: Y = L*L' ---------------
	// ********* assumed that horizon is <= 50
	compL();		// L*L' = Y
	compDnu();		// Y * dnu = -beta
	compDz();		// Phi * dz = -r_d - C'*dnu;
}

// ------- computes Phi = 2*H + kappa*P'*diag(d)^2*P ------------------
// ------------ works ------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compPhi()
{
	// ------------ first, partition diag(d) into (2N+1) diagonal matrices called d_diag[2N+1]
	// what if _nF_xTheta is not big enough?
	// Matrix<Type, Dynamic, Dynamic> d_diag[2*_N+1];
	DiagonalMatrix<Type,Dynamic> d_diag[2*_N+1];
	// int tmp = 0;
	for (int i=0 ; i <= _N-1 ; i++)
	{
		d_diag[2*i].diagonal() = d.template segment<_nInp>(i*(_nInp + _nSt));
		// tmp = tmp + Fu_rows[i];
		d_diag[2*i+1].diagonal() = d.template segment<_nSt>(_nInp+i*(_nInp + _nSt));
		// tmp = tmp + Fx_rows[i];
	}
	d_diag[2*_N].diagonal() = d.template segment<_nF_xTheta>(_N*(_nInp + _nSt));
	
	// --------------- compute elements of Phi, i.e. elements of Omicron, Pi = 2*Q_tilde, Rho, Sigma	
	
	// DiagonalMatrix<Type, _m> eyeM;
	Matrix<Type, _m, _m> eyeM;
	eyeM.setIdentity();
	// DiagonalMatrix<Type, _n> eyeN;
	Matrix<Type, _n, _n> eyeN;
	eyeN.setIdentity();
	
	// special treatment at the beginning
	// ***** note: if d_diag * d_diag is expensive, then this can also be done by building d_diag2 (squared) by computing d_i and d_i*d_i in compD()
	Omicron[0] = 2*R + kappa*Fu_transp[0]*d_diag[0]*d_diag[0]*Fu[0]	+ reg * eyeM;
	
	// cout << endl << "------- special ------" << endl;
	// cout << "reg:" << reg << endl << endl;
	// cout << "2*R + kappa*Fu_transp[0]*d_diag[0]*d_diag[0]*Fu[0]" << endl << 2*R + kappa*Fu_transp[0]*d_diag[0]*d_diag[0]*Fu[0] << endl << endl;
	// cout << "2*R + kappa*Fu_transp[0]*d_diag[0]*d_diag[0]*Fu[0] + reg * eyeM" << endl << 2*R + kappa*Fu_transp[0]*d_diag[0]*d_diag[0]*Fu[0] + reg * eyeM << endl << endl;
	// cout << "-----------" << endl << endl;
	
	// do the rest by computing three block and three block
	// might also integrate the computation of d_diag directly into this loop
	for (int i=1; i <= _N-1; i++)
	{
		if (i != _pos_omega)
		{
			Rho[i-1] = 2*Q_bar + kappa * ( Fx_transp[i-1]*d_diag[2*(i-1)+1]*d_diag[2*(i-1)+1]*Fx[i-1]  + Fu_bar_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu_bar[i] ) 
				+ reg * eyeN;
		}
		else	// i == _pos_omega
		{
			Rho[_pos_omega-1] = 2*Q_bar + kappa * ( Fx_transp[_pos_omega-1]*d_diag[2*(_pos_omega-1)+1]*d_diag[2*(_pos_omega-1)+1]*Fx[i-1]  + Fu_bar_transp[_pos_omega]*d_diag[2*_pos_omega]*d_diag[2*_pos_omega]*Fu_bar[_pos_omega] 
				+ F_xTheta_transp * d_diag[2*_N]*d_diag[2*_N]*F_xTheta)
				+ reg * eyeN;
		}
		Sigma[i-1] = 2*S + kappa * ( Fu_bar_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu[i] );
		Omicron[i] = 2*R + kappa * ( Fu_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu[i] ) 
				+ reg * eyeM;
	}
	
	// special treatment for last block
	if (_pos_omega == _N)
	{
		Rho[_N-1] = kappa * ( Fx_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx[_N-1] + F_xTheta_transp*d_diag[2*_N]*d_diag[2*_N]*F_xTheta )
				 + reg * eyeN;
	}
	else	// considered in loop above
	{
		Rho[_N-1] = kappa * ( Fx_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx[_N-1])
			 + reg * eyeN;
	}
	Sigma[_N-1] = kappa * ( F_xTheta_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta );	// independent of _pos_omega, represents the off-diag matrix
	Omicron[_N] = kappa * ( F_theta_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta )
		 + reg * eyeM;
	
	/*
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
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compPhi_tilde()
{
	Matrix<Type, _n, _n> eyeN;
	eyeN.setIdentity();
	
	// decompose Phi = L*L'
	LOmicron_diag[0].compute(Omicron[0]);
	LOmicron_diag_transp[0] = LOmicron_diag[0].matrixLLT().transpose();
	for (int i=1 ; i<= _N-1; i++)
	{
		LPi_diag[i-1].compute(2*Q_tilde);
		LPi_diag_transp[i-1] = LPi_diag[i-1].matrixLLT().transpose();
		LRho_diag[i-1].compute(Rho[i-1]);
		LRho_diag_transp[i-1] = LRho_diag[i-1].matrixLLT().transpose();
		LSigma_offDiag_transp[i-1] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(Sigma[i-1]);
		LSigma_offDiag[i-1] = LSigma_offDiag_transp[i-1].transpose();
		LOmicron_diag[i].compute(Omicron[i]-LSigma_offDiag[i-1]*LSigma_offDiag_transp[i-1]);
		LOmicron_diag_transp[i] = LOmicron_diag[i].matrixLLT().transpose();
		
		if(i == _pos_omega)
		{
			LLambda0_transp = LRho_diag[_pos_omega-1].matrixLLT().triangularView<Lower>().solve(Sigma[_N-1]);
			LLambda0 = LLambda0_transp.transpose();
			LLambda1_transp = LOmicron_diag[_pos_omega].matrixLLT().triangularView<Lower>().solve(Matrix<Type,_m,_m>::Zero() - LSigma_offDiag[_pos_omega-1]*LLambda0_transp);
			LLambda1 = LLambda1_transp.transpose();
		}
	}

	LPi_diag[_N-1].compute(2*Q_tilde_f);
	LPi_diag_transp[_N-1] = LPi_diag[_N-1].matrixLLT().transpose();
	LRho_diag[_N-1].compute(Rho[_N-1]);
	LRho_diag_transp[_N-1] = LRho_diag[_N-1].matrixLLT().transpose();

	if (_N == _pos_omega)
	{
		LSigma_offDiag_transp[_N-1] = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(Sigma[_N-1]);
		LSigma_offDiag[_N-1] = LSigma_offDiag_transp[_N-1].transpose();
		LOmicron_diag[_N].compute(Omicron[_N]-LSigma_offDiag[_N-1]*LSigma_offDiag_transp[_N-1]);
		LOmicron_diag_transp[_N] = LOmicron_diag[_N].matrixLLT().transpose();
	}
	else	// LSigma_offDiag[_N-1] is not used
	{	
		LOmicron_diag[_N].compute(Omicron[_N] - LLambda0*LLambda0_transp - LLambda1*LLambda1_transp);
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
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compY()
{
	// computation of Y is done in three steps: Y = C * X
	// 1. L*U = C'
	// 2. L'*X = U --> L*L'*X = Phi*X = C', i.e. X = Phi_tilde * C'
	// 3. C*X = C*Phi_tilde*C' = Y

	Matrix<Type, _n, _n> eye;
	eye.setIdentity();
	// Matrix<Type, _m, _n> zero;
	// zero.setZero();

	// 1. Compute elements of Matrix U
	// treat the first 2 U_bar and first 3 U specially
	U[0] = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(-Bm_transp);
	U_bar[0] = LPi_diag[0].matrixLLT().triangularView<Lower>().solve(eye);
	
	U[1] = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(-B_transp);
	U_bar[1] = LRho_diag[0].matrixLLT().triangularView<Lower>().solve(eye);
	U[2] = LOmicron_diag[1].matrixLLT().triangularView<Lower>().solve( /*zero*/ - LSigma_offDiag[0]*U_bar[1] );
	
	// remaining U_bar and U in middle have structures
	for (int i=1; i<= _N-2; i++)
	{
		U_bar[2+(i-1)*5] = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(-Am_tilde_transp);
		U_bar[3+(i-1)*5] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(-Bm_bar_transp);
		U[3+(i-1)*3] = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(-Bm_transp - LSigma_offDiag[i-1]*U_bar[3+(i-1)*5] );
		U_bar[4+(i-1)*5] = LPi_diag[i].matrixLLT().triangularView<Lower>().solve( eye );
		
		U_bar[5+(i-1)*5] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(-A_bar_transp);
		U[4+(i-1)*3] = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(-B_transp - LSigma_offDiag[i-1]*U_bar[5+(i-1)*5]);
		U_bar[6+(i-1)*5] = LRho_diag[i].matrixLLT().triangularView<Lower>().solve( eye );
		U[5+(i-1)*3] = LOmicron_diag[i+1].matrixLLT().triangularView<Lower>().solve( /*zero*/ - LSigma_offDiag[i]*U_bar[6+(i-1)*5] );
	}
	
	U_bar[2+(_N-2)*5] = LPi_diag[_N-2].matrixLLT().triangularView<Lower>().solve(-Am_tilde_transp);
	U_bar[3+(_N-2)*5] = LRho_diag[_N-2].matrixLLT().triangularView<Lower>().solve(-Bm_bar_transp);
	U[3+(_N-2)*3] = LOmicron_diag[_N-1].matrixLLT().triangularView<Lower>().solve(-Bm_transp - LSigma_offDiag[_N-2]*U_bar[3+(_N-2)*5] );
	U_bar[4+(_N-2)*5] = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve( eye );
	
	U_bar[5+(_N-2)*5] = LRho_diag[_N-2].matrixLLT().triangularView<Lower>().solve(-A_bar_transp);
	U[4+(_N-2)*3] = LOmicron_diag[_N-1].matrixLLT().triangularView<Lower>().solve(-B_transp - LSigma_offDiag[_N-2]*U_bar[5+(_N-2)*5]);
	U_bar[6+(_N-2)*5] = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve( eye );
	
	if (_N == _pos_omega)
	{
		U[5+(_N-2)*3] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( /*zero*/ - LSigma_offDiag[_N-1]*U_bar[6+(_N-2)*5] );
	}
	else
	{
		UO[0] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(/*zero*/ - LLambda0*U_bar[1+(_pos_omega-1)*5] - LLambda1*U[2+(_pos_omega-1)*3]);
		UO[1] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(/*zero*/ - LLambda0*U_bar[3+(_pos_omega-1)*5] - LLambda1*U[3+(_pos_omega-1)*3]);
		UO[2] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(/*zero*/ - LLambda0*U_bar[5+(_pos_omega-1)*5] - LLambda1*U[4+(_pos_omega-1)*3]);
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
		XO[0] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[0]);
		XO[1] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[1]);
		XO[2] = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(UO[2]);
	}
	
	X[0] = LOmicron_diag_transp[0].template triangularView<Upper>().solve(U[0]);
	X_bar[0] = LPi_diag_transp[0].template triangularView<Upper>().solve(U_bar[0]);
	
	if (_pos_omega == 1)	// build col1
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
		if (i == _pos_omega)	// col2
		{
			X_bar[2+(_pos_omega-1)*5] = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[2+(_pos_omega-1)*5]);
			X[3+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(U[3+(_pos_omega-1)*3]-LLambda1_transp*XO[1]);
			X_bar[3+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[3+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[3+(_pos_omega-1)*3] - LLambda0_transp*XO[1]);
			X_bar[4+(_pos_omega-1)*5] = LPi_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar[4+(_pos_omega-1)*5]);
		}
		else	// standard
		{
			X_bar[2+(i-1)*5] = LPi_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[2+(i-1)*5]);
			X[3+(i-1)*3] = LOmicron_diag_transp[i].template triangularView<Upper>().solve(U[3+(i-1)*3]);
			X_bar[3+(i-1)*5] = LRho_diag_transp[i-1].template triangularView<Upper>().solve(U_bar[3+(i-1)*5] - LSigma_offDiag_transp[i-1]*X[3+(i-1)*3]);
			X_bar[4+(i-1)*5] = LPi_diag_transp[i].template triangularView<Upper>().solve(U_bar[4+(i-1)*5]);
		}
		
		if (i == _pos_omega)	// col3
		{
			X[4+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve( U[4+(_pos_omega-1)*3] - LLambda1_transp*XO[2]);
			X_bar[5+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[5+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[4+(_pos_omega-1)*3] - LLambda0_transp*XO[2]);
			X[5+(_pos_omega-1)*3] = LOmicron_diag_transp[_pos_omega+1].template triangularView<Upper>().solve(U[5+(_pos_omega-1)*3]);
			X_bar[6+(_pos_omega-1)*5] = LRho_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar[6+(_pos_omega-1)*5] - LSigma_offDiag_transp[_pos_omega]*X[5+(_pos_omega-1)*3]);
		}
		else if(i == _pos_omega-1)	// col1
		{
			X[4+(_pos_omega-2)*3] = LOmicron_diag_transp[_pos_omega-1].template triangularView<Upper>().solve( U[4+(_pos_omega-2)*3] );
			X_bar[5+(_pos_omega-2)*5] = LRho_diag_transp[_pos_omega-2].template triangularView<Upper>().solve(U_bar[5+(_pos_omega-2)*5] - LSigma_offDiag_transp[_pos_omega-2]*X[4+(_pos_omega-2)*3]);
			X[5+(_pos_omega-2)*3] = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(U[5+(_pos_omega-2)*3] - LLambda1_transp*XO[0]);
			X_bar[6+(_pos_omega-2)*5] = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar[6+(_pos_omega-2)*5] - LSigma_offDiag_transp[_pos_omega-1]*X[5+(_pos_omega-2)*3] - LLambda0_transp*XO[0]);		
		}
		else	// if the off-diag element in P*z<h has no influence
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
	else if(_pos_omega == _N-1)	// compute col2 and col3
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


// ------------ function computes beta = -r_p + C*Phi_tilde*r_d --------
// ------------ works --------------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compBeta()
{
	// solve in steps
	// beta = -r_p + C*tmp1, where tmp1: Phi*tmp1 = r_d
	// 1. L*tmp2 = r_d   --> compute tmp2
	// 2. L'*tmp1 = tmp2  ---> compute tmp1
	// 3. beta = -r_p + C*tmp1
	
	
	// 1. compute tmp2: L*tmp2 = r_d
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> tmp2;
	tmp2.template segment<_m>(0) = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_m>(0) );
	for (int i = 1; i <= _N-1; i++)
	{
		tmp2.template segment<_n>(_m+(i-1)*offset) = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_n>(_m+(i-1)*offset));
		tmp2.template segment<_n>(_m+(i-1)*offset+_n) = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_n>(_m+(i-1)*offset+_n));
		tmp2.template segment<_m>(_m+(i-1)*offset+_n+_n) = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_m>(_m+(i-1)*offset+_n+_n) - LSigma_offDiag[i-1]*tmp2.template segment<_n>(_m+(i-1)*offset+_n));
	}
	if (_pos_omega == _N)
	{
		tmp2.template segment<_n>(_m+(_N-1)*offset) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_n>(_m+(_N-1)*offset));
		tmp2.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_n>(_m+(_N-1)*offset+_n));
		tmp2.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n) - LSigma_offDiag[_N-1]*tmp2.template segment<_n>(_m+(_N-1)*offset+_n));
	}
	else	// need to use LLambda0 and LLambda1
	{
		tmp2.template segment<_n>(_m+(_N-1)*offset) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_n>(_m+(_N-1)*offset));
		tmp2.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_n>(_m+(_N-1)*offset+_n));
		tmp2.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n) - LLambda0*tmp2.template segment<_n>(_m+(_pos_omega-1)*offset+_n) - LLambda1*tmp2.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) );	
	}
	/*
	cout << "tmp2:" << setprecision(15)<< endl << tmp2 << endl << endl;
	cout << "r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n)" << endl << r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n) << endl << endl;
	cout << "LOmicron_diag[_N]" << endl << LOmicron_diag[_N].matrixLLT() << endl << endl;
	cout << "LLambda0" << endl << LLambda0 << endl << endl;
	cout << "tmp2.template segment<_n>(_m+(_pos_omega-1)*offset+_n)" << endl << tmp2.template segment<_n>(_m+(_pos_omega-1)*offset+_n) << endl << endl;
	cout << "LLambda1" << endl << LLambda1 << endl << endl;
	cout << "tmp2.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n)" << endl << tmp2.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) << endl << endl;
	*/
	
	
	// 2. compute tmp1: L'*tmp1 = tmp2
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> tmp1;
	tmp1.template segment<_m>(0) = LOmicron_diag_transp[0].template triangularView<Upper>().solve(tmp2.template segment<_m>(0));
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
	else	// standard and compute missing block
	{
		tmp1.template segment<_n>(_m+(_N-1)*offset) = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_N-1)*offset));
		tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(_N-1)*offset+_n+_n));
		tmp1.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_N-1)*offset+_n));
		
		tmp1.template segment<_n>(_m+(_pos_omega-1)*offset) = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_pos_omega-1)*offset));
		tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp2.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda1_transp*tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) );
		tmp1.template segment<_n>(_m+(_pos_omega-1)*offset+_n) = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2.template segment<_n>(_m+(_pos_omega-1)*offset+_n) - LSigma_offDiag_transp[_pos_omega-1]*tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda0_transp*tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) );
	}
		
	// cout << "tmp1:" << endl << tmp1 << endl << endl;
	
	// 3. beta = -r_p + C*tmp1
	beta.template segment<_n>(0) = -r_p.template segment<_n>(0) + ( -Bm*tmp1.template segment<_m>(0) + tmp1.template segment<_n>(_m) );
	beta.template segment<_n>(_n) = -r_p.template segment<_n>(_n) + (- B*tmp1.template segment<_m>(0) + tmp1.template segment<_n>(_m+_n) );

	for(int i=1; i<= _N-1; i++)
	{
		beta.template segment<_n>(2*i*_n) = -r_p.template segment<_n>(2*i*_n) + (    
			- Am_tilde*tmp1.template segment<_n>(_m+(i-1)*offset) 
			- Bm_bar * tmp1.template segment<_n>(_m+(i-1)*offset+_n)
			- Bm* tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n) 
			+ tmp1.template segment<_n>(i*offset+_m)    );
					
		beta.template segment<_n>(2*i*_n+_n) = -r_p.template segment<_n>( 2*i*_n+_n) + (  
		    - A_bar * tmp1.template segment<_n>(_m+_n+(i-1)*offset)
		    - B * tmp1.template segment<_m>(_m+_n+(i-1)*offset+_n)
			+ tmp1.template segment<_n>(_m+(i)*offset+_n)        );
	}
	// cout << "beta:" << endl << beta << endl << endl;
}


// ------------ function computes L: L*L' = Y --------------------------
// -------------- compute components of L[3][2*_N] ---------------------
// ---- remark: L[i][i] are lower triangular matrices and can be accessed by L[i][i].matrixLLT().triangularView<Lower>()
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compL()		// L*L' = Y
{
	// special treatment for L11, L21, L31
	L_diag[0].compute(Y[0][0]);		// Y11 = L11*L11'
	L_diag_transp[0] = L_diag[0].matrixLLT().transpose();
	
	L_offDiag_transp[0][0] = L_diag[0].matrixLLT().triangularView<Lower>().solve(Y[0][1]);
	L_offDiag[0][0] = L_offDiag_transp[0][0].transpose();
	
	L_offDiag_transp[1][0] = L_diag[0].matrixLLT().triangularView<Lower>().solve(Y[0][2]);
	L_offDiag[1][0] = L_offDiag_transp[1][0].transpose();
	//cout <<  "L_diag[0]:" << endl << L_diag[0].matrixLLT() << endl << endl;
	//cout << "L_offDiag[0][0]:" << endl << L_offDiag[0][0] << endl << endl;
	//cout << "L_offDiag[1][0]:" << endl << L_offDiag[1][0] << endl << endl;
	
	// special treatment for L22, L32, L42
	L_diag[1].compute( Y[1][1]-L_offDiag[0][0]*L_offDiag_transp[0][0] );
	L_diag_transp[1] = L_diag[1].matrixLLT().transpose();
	
	L_offDiag_transp[0][1] = L_diag[1].matrixLLT().triangularView<Lower>().solve(Y[1][2]-L_offDiag[0][0]*L_offDiag_transp[1][0]);
	L_offDiag[0][1] = L_offDiag_transp[0][1].transpose();
	
	L_offDiag_transp[1][1] = L_diag[1].matrixLLT().triangularView<Lower>().solve(Y[0][3]);
	L_offDiag[1][1] = L_offDiag_transp[1][1].transpose();
	// cout <<  "L_diag[1]:" << endl << L_diag[1].matrixLLT() << endl << endl;
	// cout << "L_offDiag[0][1]:" << endl << L_offDiag[0][1] << endl << endl;
	// cout << "L_offDiag[1][1]:" << endl << L_offDiag[1][1] << endl << endl;
	
	// cases in the middle
	for (int i = 1; i <= 2*_N-4; i++)
	{
		L_diag[i+1].compute( Y[2][i+1] - L_offDiag[1][i-1]*L_offDiag_transp[1][i-1] - L_offDiag[0][i]*L_offDiag_transp[0][i] );
		L_diag_transp[i+1] = L_diag[i+1].matrixLLT().transpose();
		L_offDiag_transp[0][i+1] = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( Y[1][i+2] - L_offDiag[0][i]*L_offDiag_transp[1][i] );
		L_offDiag[0][i+1] = L_offDiag_transp[0][i+1].transpose();
		
		L_offDiag_transp[1][i+1] = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( Y[0][i+3] );
		L_offDiag[1][i+1] = L_offDiag_transp[1][i+1].transpose();
	//	 cout <<  "L_diag[i+1]:" << endl << L_diag[i+1].matrixLLT() << endl << endl;
	//	 cout << "L_offDiag[0][i+1]:" << endl << L_offDiag[0][i+1] << endl << endl;
	//	 cout << "L_offDiag[1][i+1]:" << endl << L_offDiag[1][i+1] << endl << endl;
	}
		
	
	// special treatment in the end, i.e. i = 2*_N-3
	L_diag[2*_N-2].compute( Y[2][2*_N-2] - L_offDiag[1][2*_N-4]*L_offDiag_transp[1][2*_N-4] - L_offDiag[0][2*_N-3]*L_offDiag_transp[0][2*_N-3] );
	L_diag_transp[2*_N-2] = L_diag[2*_N-2].matrixLLT().transpose();
	L_offDiag_transp[0][2*_N-2] = L_diag[2*_N-2].matrixLLT().triangularView<Lower>().solve( Y[1][2*_N-1] - L_offDiag[0][2*_N-3]*L_offDiag_transp[1][2*_N-3] );
	L_offDiag[0][2*_N-2] = L_offDiag_transp[0][2*_N-2].transpose();
	
	// i = 2*_N-2
	L_diag[2*_N-1].compute( Y[2][2*_N-1] - L_offDiag[1][2*_N-3]*L_offDiag_transp[1][2*_N-3] - L_offDiag[0][2*_N-2]*L_offDiag_transp[0][2*_N-2] );
	L_diag_transp[2*_N-1] = L_diag[2*_N-1].matrixLLT().transpose();	
}


// ------------ function computes L*L'*dnu = -beta ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDnu()
{
	// 1) first, solve for delta: L*delta = -beta
	Matrix<Type, 2*_N*_n, 1> delta;
	
	// special cases in the beginning
	delta.template segment<_n>(0) = L_diag[0].matrixLLT().triangularView<Lower>().solve(-beta.template segment<_n>(0));
	delta.template segment<_n>(_n) = L_diag[1].matrixLLT().triangularView<Lower>().solve(-beta.template segment<_n>(_n) - L_offDiag[0][0]*delta.template segment<_n>(0));
	
		
	// remaining cases are regular
	for (int i=1; i<= 2*_N-2; i++)
	{
		// delta.segment(_n+i*_n,_n) = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( -beta.segment(_n+i*_n,_n) - L_offDiag[1][i-1]*delta.segment((i-1)*_n,_n) - L_offDiag[0][i]*delta.segment(i*_n,_n) );
		delta.template segment<_n>(_n+i*_n) = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( -beta.template segment<_n>(_n+i*_n) - L_offDiag[1][i-1]*delta.template segment<_n>((i-1)*_n) - L_offDiag[0][i]*delta.template segment<_n>(i*_n) );
	}
	//cout << "delta:" << endl << delta << endl << endl;
	
	
	// 2) now, solve for L'*Dnu = delta
	dnu.template segment<_n>(2*_n*_N - _n) = L_diag_transp[2*_N-1].template triangularView<Upper>().solve(delta.template segment<_n>(2*_n*_N - _n) );
	dnu.template segment<_n>(2*_n*_N - _n - _n) = L_diag_transp[2*_N-2].template triangularView<Upper>().solve( delta.template segment<_n>(2*_n*_N - _n - _n) - L_offDiag_transp[0][2*_N-2]*dnu.template segment<_n>(2*_n*_N - _n) );
	
	//remaining cases are regular
	for (int i=1; i<=2*_N-2; i++)
	{
		// dnu.segment(2*_n*_N-(i+2)*_n,_n) = L_diag[2*_N-(i+2)].matrixLLT().transpose().triangularView<Upper>().solve( delta.segment(2*_n*_N-(i+2)*_n,_n) - L_offDiag_transp[0][2*_N-(i+2)]*dnu.segment(2*_n*_N-(i+1)*_n,_n) - L_offDiag_transp[1][2*_N-(i+2)]*dnu.segment(2*_n*_N-i*_n,_n)  );
		dnu.template segment<_n>(2*_n*_N-(i+2)*_n) = L_diag_transp[2*_N-(i+2)].template triangularView<Upper>().solve( delta.template segment<_n>(2*_n*_N-(i+2)*_n) - L_offDiag_transp[0][2*_N-(i+2)]*dnu.template segment<_n>(2*_n*_N-(i+1)*_n) - L_offDiag_transp[1][2*_N-(i+2)]*dnu.template segment<_n>(2*_n*_N-i*_n)  );
	}
	//cout << "dnu" << endl << dnu << endl << endl;
}


// ------------ function computes Phi * dz = -r_d - C' * dnu ---------------------
// ------------ implementing --------------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDz()
{	
	// computed in two parts
	// 1. tmp = -r_d - C' * dnu
	// 2. L*L'*dz = tmp
	// 3. L*tmp1 = tmp
	// 4. L'*dz = tmp1
	
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> tmp;
	tmp.template segment<_m>(0) = -r_d.template segment<_m>(0) + Bm_transp*dnu.template segment<_n>(0) + B_transp*dnu.template segment<_n>(_n);
	
	for (int i=1; i<= _N-1; i++)
	{
		tmp.template segment<_n>(_m+(i-1)*offset) = -r_d.template segment<_n>(_m+(i-1)*offset) - dnu.template segment<_n>(2*(i-1)*_n) + Am_tilde_transp*dnu.template segment<_n>(2*(i-1)*_n+_n+_n) ;
		tmp.template segment<_n>(_m+(i-1)*offset+_n) = -r_d.template segment<_n>(_m+(i-1)*offset+_n) - dnu.template segment<_n>(2*(i-1)*_n+_n) + Bm_bar_transp*dnu.template segment<_n>(2*(i-1)*_n+_n+_n) + A_bar_transp * dnu.template segment<_n>(2*(i-1)*_n+_n+_n+_n);
		tmp.template segment<_m>(_m+(i-1)*offset+_n+_n) = -r_d.template segment<_m>(_m+(i-1)*offset+_n+_n) + Bm_transp*dnu.template segment<_n>(2*(i-1)*_n+_n+_n) + B_transp*dnu.template segment<_n>(2*(i-1)*_n+_n+_n+_n);
	}
	
	tmp.template segment<_n>(_m+(_N-1)*offset) = -r_d.template segment<_n>(_m+(_N-1)*offset) - dnu.template segment<_n>(2*(_N-1)*_n);
	tmp.template segment<_n>(_m+(_N-1)*offset+_n) = -r_d.template segment<_n>(_m+(_N-1)*offset+_n) - dnu.template segment<_n>(2*(_N-1)*_n+_n);
	tmp.template segment<_m>(_m+(_N-1)*offset+_n+_n) = -r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n);
	
	//cout << "tmp:" << endl << tmp << endl << endl;
	
	
	// 3. L*tmp1 = tmp
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> tmp1;
	tmp1.template segment<_m>(0) = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_m>(0) );
	for (int i = 1; i <= _N-1; i++)
	{
		tmp1.template segment<_n>(_m+(i-1)*offset) = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_n>(_m+(i-1)*offset));
		tmp1.template segment<_n>(_m+(i-1)*offset+_n) = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_n>(_m+(i-1)*offset+_n));
		tmp1.template segment<_m>(_m+(i-1)*offset+_n+_n) = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_m>(_m+(i-1)*offset+_n+_n) - LSigma_offDiag[i-1]*tmp1.template segment<_n>(_m+(i-1)*offset+_n));
	}
	if (_pos_omega == _N)
	{
		tmp1.template segment<_n>(_m+(_N-1)*offset) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_n>(_m+(_N-1)*offset));
		tmp1.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_n>(_m+(_N-1)*offset+_n));
		tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_m>(_m+(_N-1)*offset+_n+_n) - LSigma_offDiag[_N-1]*tmp1.template segment<_n>(_m+(_N-1)*offset+_n));
	}
	else
	{
		tmp1.template segment<_n>(_m+(_N-1)*offset) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_n>(_m+(_N-1)*offset));
		tmp1.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_n>(_m+(_N-1)*offset+_n));
		tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(tmp.template segment<_m>(_m+(_N-1)*offset+_n+_n) - LLambda0*tmp1.template segment<_n>(_m+(_pos_omega-1)*offset+_n) - LLambda1*tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n)  );
	}
	// cout << "tmp1:" << endl << tmp1 << endl << endl;
	
	
	// 4. L'*dz = tmp1
	dz.template segment<_m>(0) = LOmicron_diag_transp[0].template triangularView<Upper>().solve(tmp1.template segment<_m>(0));
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
	else	// standard ending and missing block
	{
		dz.template segment<_n>(_m+(_N-1)*offset) = LPi_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_N-1)*offset));
		dz.template segment<_m>(_m+(_N-1)*offset+_n+_n) = LOmicron_diag_transp[_N].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(_N-1)*offset+_n+_n));
		dz.template segment<_n>(_m+(_N-1)*offset+_n) = LRho_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_N-1)*offset+_n) );
	
		dz.template segment<_n>(_m+(_pos_omega-1)*offset) = LPi_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_pos_omega-1)*offset));
		dz.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) = LOmicron_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp1.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda1_transp*dz.template segment<_m>(_m+(_N-1)*offset+_n+_n)  );
		dz.template segment<_n>(_m+(_pos_omega-1)*offset+_n) = LRho_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1.template segment<_n>(_m+(_pos_omega-1)*offset+_n) - LSigma_offDiag_transp[_pos_omega-1]*dz.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) - LLambda0_transp*dz.template segment<_m>(_m+(_N-1)*offset+_n+_n)  );
	}
	// cout << "dz" << endl << dz << endl << endl;
}


// ------------ function tests if z_warm satisfies P*z < h --------
// ------- might consider computing entire matrix check and then iterate for neg. elements ----
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
double LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compT()
{
	Matrix<Type, _m, 1> c_tmp;		// use this to temporarily copy c-blocks out of z
	Matrix<Type, _m, 1> dc_tmp;		// use this to temporarily copy c-blocks out of dz
	Matrix<Type, _n, 1> x_bar_tmp;	// use this to temporarily copy x_bar-blocks out of z
	Matrix<Type, _n, 1> dx_bar_tmp;	// use this to temporarily copy x_bar-blocks out of dz
	Matrix<Type, Dynamic, 1> check;		// length of vector not yet fixed
	Matrix<Type, Dynamic, 1> dcheck;	// dcheck = P*dz
	Matrix<Type, _nSt, 1> t_vec_nSt;		// t_vec = check ./ dcheck
	Matrix<Type, _nInp, 1> t_vec_nInp;
	Matrix<Type, _nF_xTheta, 1> t_vec_nF_xTheta;
	
	double t = 1;	// stores smallest variable
	double t_tmp;	// t_tmp = min(t_vec)
	
	// cout << endl << "======= starting testPos() =======" << endl;
	// special treatment at beginning
	c_tmp = z.template segment<_m>(0);
	dc_tmp = dz.template segment<_m>(0);
	check = (fu[0] - Fu_bar[0]*x_hat) - (Fu[0]*c_tmp);	// should be >0
	dcheck = Fu[0] * dc_tmp;
	// t_vec.resize(_nInp);
	t_vec_nInp.setConstant(1);
	for (int j=1; j <= _nInp; j++)
	{
		if (dcheck[j-1] > 0)	// neg. cases not interesting
		{
			t_vec_nInp[j-1] = check[j-1]/dcheck[j-1];
		}
	}
	t_tmp = t_vec_nInp.minCoeff();
	if (t_tmp < t)
	{
		t = t_tmp;
	}

	// cout << "finished special treatment at start" << endl;
	
	// class variable: offset = _n + _n + _m;
	// general treatment in the middle 
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	// compute (h - P*z)_i
		// cout << "round " << i << "in general treatment" << endl;
		x_bar_tmp = z.template segment<_n>((i-1)*offset+_m+_n);
		dx_bar_tmp = dz.template segment<_n>((i-1)*offset+_m+_n);
		c_tmp = z.template segment<_m>((i-1)*offset+_m+_n+_n);
		dc_tmp = dz.template segment<_m>((i-1)*offset+_m+_n+_n);
		check = fx[i-1] - Fx[i-1]*x_bar_tmp;
		dcheck = Fx[i-1]*dx_bar_tmp;
		// t_vec.resize(_nSt);
		t_vec_nSt.setConstant(1);
		// cout << "computed check" << endl << check << endl << endl;
		for (int j=1; j<= _nSt; j++)
		{
			if(dcheck[j-1]>0)
			{
				t_vec_nSt[j-1] = check[j-1]/dcheck[j-1];
			}
		}
		t_tmp = t_vec_nSt.minCoeff();
		if (t_tmp < t)
		{
			t = t_tmp;
		}
		
		check = fu[i] - (Fu_bar[i]*x_bar_tmp + Fu[i]*c_tmp);
		dcheck = Fu_bar[i]*dx_bar_tmp + Fu[i]*dc_tmp;
		// t_vec.resize(_nInp);
		t_vec_nInp.setConstant(1);
		for (int j=1; j<=_nInp; j++)
		{
			if(dcheck[j-1]>0)
			{
				t_vec_nInp[j-1] = check[j-1]/dcheck[j-1];
			}
		}
		t_tmp = t_vec_nInp.minCoeff();
		if (t_tmp < t)
		{
			t = t_tmp;
		}
		//cout << "for-loop completion of round" << i << endl;
	}
	
	// cout << "finished general treatment" << endl;
	
	// special case for last blocks
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);
	dx_bar_tmp = dz.template segment<_n>((_N-1)*offset+_m+_n);
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	dc_tmp = dz.template segment<_m>((_N-1)*offset+_m+_n+_n);
	check = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	dcheck = Fx[_N-1]*dx_bar_tmp;
	// t_vec.resize(_nSt);
	t_vec_nSt.setConstant(1);
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j<= _nSt; j++)
	{
		if (dcheck[j-1]>0)
		{
			t_vec_nSt[j-1] = check[j-1]/dcheck[j-1];
		}
	}
	t_tmp = t_vec_nSt.minCoeff();
	if (t_tmp < t)
	{
		t = t_tmp;
	}
	
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);
	dx_bar_tmp = dz.template segment<_n>((_pos_omega-1)*offset+_m+_n);
	check = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	dcheck = F_xTheta*dx_bar_tmp + F_theta*dc_tmp;
	// t_vec_nF_xTheta.resize(_nF_xTheta);
	t_vec_nF_xTheta.setConstant(1);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if(dcheck[j-1]>0)
		{
			t_vec_nF_xTheta[j-1] = check[j-1]/dcheck[j-1];
		}
	}
	t_tmp = t_vec_nF_xTheta.minCoeff();
	if (t_tmp < t)
	{
		t = t_tmp;
	}
	
	// return the result
	if (t == 1)
	{
		return 1;
	} 
	else
	{
		return 0.99*t;	// guarantees strict feasibility
	}
}


// ------------ function computes cost vectors q_tilde_vec, q_bar_vec, r_vec --------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compRQ()
{
	// compute the vectors u_star
	// "ColPivHouseholderQRPreconditioner" is more accurate, "HouseholderQRPreconditioner" is faster
	JacobiSVD<MatrixXd, HouseholderQRPreconditioner> svd(Bm, ComputeThinU | ComputeThinV);	// computes SVD of Bm
	u_star[0] = svd.solve(x_star[0] - Am_tilde*x_hat - tm_tilde);	
	for (int i=1; i <= _N-1; i++)
	{
		u_star[i] = svd.solve(x_star[i] - Am_tilde*x_star[i-1] - tm_tilde);
	}
	
	/*
	for (int i = 0; i <= _N-1; i++)
	{
		cout << "u_star[" << i << "]" << endl << u_star[i] << endl << endl;
	}
	*/
	
	// compute the vectors q_tilde[]; q_tilde_f; r[]
	for (int i = 0; i <= _N-1; i++)
	{
		q_tilde_vec[i] = -2*x_star[i];
		// q_tilde_vec[i] << 1 , 1 , 2 , -4 , 5;
		r_vec[i] = -2*u_star[i];
		// r_vec[i] << 1, 1;
		q_bar_vec[i] = K_transp*r_vec[i];
	}
	
	/*
	for (int i = 0; i <= _N-1; i++)
	{
		cout << "q_tilde_vec[" << i << "]" << endl << q_tilde_vec[i] << endl << endl;
		cout << "r_vec[" << i << "]" << endl << r_vec[i] << endl << endl;
		cout << "q_bar_vec[" << i << "]" << endl << q_bar_vec[i] << endl << endl;
	}
	*/
	
}



// ------------ function computes augmented matrices needed by phaseI ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compMatrices_PhaseI()
{
	// vvv build the matrices Fx_hat(_transp), Fu_hat(_transp), F_xTheta_hat(_transp), F_theta_hat(_transp), Fu_bar(_transp)
	Matrix<Type, _nInp, _nInp> eye_nInp;
	eye_nInp.setIdentity();
	Matrix<Type, _nSt, _nSt> eye_nSt;
	eye_nSt.setIdentity();
	double weight = 100;
	
	// Fu_hat[0].resize(_nInp,_m+_nInp);
	// Fu_hat[0] << -Fu[0], Matrix<Type, Dynamic, Dynamic, 0, _nInp, _nInp>::Identity(_nInp,_nInp);
	Fu_hat[0] << -Fu[0], eye_nInp;
	Fu_hat[0] = -Fu_hat[0];
	Fu_hat_transp[0] = Fu_hat[0].transpose();
	for (int i=1; i<=_N-1; i++)
	{
		// Fx_hat[i-1].resize(_nSt,_n+_nSt);
		// Fx_hat[i-1] << -Fx[i-1], Matrix<Type, Dynamic, Dynamic, 0, _nSt, _nSt>::Identity(_nSt,_nSt);
		Fx_hat[i-1] << -Fx[i-1], eye_nSt;
		Fx_hat[i-1] = -Fx_hat[i-1];
		Fx_hat_transp[i-1] = Fx_hat[i-1].transpose();
		
		// Fu_bar_hat[i-1] = Matrix<Type, Dynamic, Dynamic, 0, _nInp, _n+_nSt>::Zero(_nInp,_n+_nSt);
		Fu_bar_hat[i-1].setZero();
		Fu_bar_hat[i-1].block(0,0,_nInp,_n) = Fu_bar[i];
		Fu_bar_hat_transp[i-1] = Fu_bar_hat[i-1].transpose();
		
		// Fu_hat[i].resize(_nInp,_m+_nInp);
		Fu_hat[i] << -Fu[i], eye_nInp;
		Fu_hat[i] = -Fu_hat[i];
		Fu_hat_transp[i] = Fu_hat[i].transpose();
	}
	// Fx_hat[_N-1].resize(_nSt,_n+_nSt);
	Fx_hat[_N-1] << -Fx[_N-1], eye_nSt;
	Fx_hat[_N-1] = -Fx_hat[_N-1];
	Fx_hat_transp[_N-1] = Fx_hat[_N-1].transpose();
	
	F_xTheta_hat.setZero();
	F_xTheta_hat.block(0,0,_nF_xTheta,_n) = F_xTheta;
	F_xTheta_hat_transp = F_xTheta_hat.transpose();
	
	F_theta_hat << -F_theta, Matrix<Type, _nF_xTheta, _nF_xTheta>::Identity(_nF_xTheta,_nF_xTheta);
	F_theta_hat = -F_theta_hat;
	F_theta_hat_transp = F_theta_hat.transpose();
	// ^^^ built the matrices
	
	// vvv definitions of Bm_hat, B_hat, B_bar_hat, A_bar_hat, Identity_hat
	for(int i=1; i<= _N-1; i++)
	{
		// Bm_hat[i-1] = Matrix<Type, _n, Dynamic, 0, _n, _m+_nInp>::Zero(_n, _m+_nInp);
		Bm_hat[i-1].setZero();
		Bm_hat[i-1].template block<_n,_m>(0,0) = Bm;
		Bm_hat_transp[i-1] = Bm_hat[i-1].transpose();
		// B_hat[i-1] = Matrix<Type, _n, Dynamic, 0, _n, _m+_nInp>::Zero(_n, _m+_nInp);
		B_hat[i-1].setZero();
		B_hat[i-1].template block<_n,_m>(0,0) = B;
		B_hat_transp[i-1] = B_hat[i-1].transpose();
		// Identity_hat[i-1] = Matrix<Type, _n, Dynamic, 0, _n, _n+_nSt>::Zero(_n, _n+_nSt);
		Identity_hat[i-1].setZero();
		Identity_hat[i-1].topLeftCorner(_n,_n).setIdentity();	// or: Identity_hat[i-1].topLeftCorner(_n,_n).setIdentity();
		Identity_hat_transp[i-1] = Identity_hat[i-1].transpose();
		// Bm_bar_hat[i-1] = Matrix<Type, _n, Dynamic, 0, _n, _n+_nSt>::Zero(_n, _n+_nSt);
		Bm_bar_hat[i-1].setZero();
		Bm_bar_hat[i-1].template block<_n,_n>(0,0) = Bm_bar;
		Bm_bar_hat_transp[i-1] = Bm_bar_hat[i-1].transpose();
		A_bar_hat[i-1].setZero();
		A_bar_hat[i-1].template block<_n,_n>(0,0) = A_bar;
		A_bar_hat_transp[i-1] = A_bar_hat[i-1].transpose();
	}
	// Bm_hat[_N-1] = Matrix<Type, _n, Dynamic, 0, _n, _m+_nInp>::Zero(_n, _m+_nInp);
	Bm_hat[_N-1].setZero();
	Bm_hat[_N-1].template block<_n,_m>(0,0) = Bm;
	Bm_hat_transp[_N-1] = Bm_hat[_N-1].transpose();
	// B_hat[_N-1] = Matrix<Type, _n, Dynamic, 0, _n, _m+_nInp>::Zero(_n, _m+_nInp);
	B_hat[_N-1].setZero();
	B_hat[_N-1].template block<_n,_m>(0,0) = B;
	B_hat_transp[_N-1] = B_hat[_N-1].transpose();
	// Identity_hat[_N-1] = Matrix<Type, _n, Dynamic, 0, _n, _n+_nSt>::Zero(_n, _n+_nSt);
	Identity_hat[_N-1].setZero();
	Identity_hat[_N-1].topLeftCorner(_n,_n).setIdentity();	// or: Identity_hat[i-1].topLeftCorner(_n,_n).setIdentity();
	Identity_hat_transp[_N-1] = Identity_hat[_N-1].transpose();
	// ^^^ ended definition of Bm_hat, B_hat, B_bar_hat, A_bar_hat, Identity_hat
	
	
	Matrix<Type, _m+_nInp,1> const_c;
	const_c.setConstant(weight);
	const_c.template segment<_m>(0).setZero();
	Matrix<Type, _n+_nSt,1> const_x;
	const_x.setConstant(weight);
	const_x.template segment<_n>(0).setZero();
	
	// vvv build elements in cost vector g_hat
	for (int i=1; i<=_N ; i++)
	{
		// g_hat_c[i-1] = Matrix<Type, Dynamic, 1, 0, _m+_nInp, 1>::Constant(_m+_nInp,100); // set to 100
		// g_hat_c[i-1].template segment<_m>(0).setZero();
		g_hat_c[i-1] = const_c;
		// g_hat_x[i-1] = Matrix<Type, Dynamic, 1, 0, _n+_nSt ,1>::Constant(_n+_nSt,100);	// set to 100
		// g_hat_x[i-1].template segment<_n>(0).setZero(); // set first couples to 0
		g_hat_x[i-1] = const_x;
	}
	// g_hat_theta = Matrix<Type, _m+_nF_xTheta, 1>::Constant(_m+_nF_xTheta,100);
	g_hat_theta.setConstant(weight);
	g_hat_theta.template segment<_m>(0).setZero();
	
	/*
	for (int i = 1; i<=_N; i++)
	{
		cout << "g_hat_c[" << i-1 << "]" << endl << g_hat_c[i-1] << endl << endl;
		cout << "g_hat_x[" << i-1 << "]" << endl << g_hat_x[i-1] << endl << endl;
	}
	cout << "g_hat_theta" << endl << g_hat_theta << endl << endl;
	*/
	
	// ^^^ elements of g_hat built
}


// ------------ function computes a feasibe gamma needed by phaseI ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compGamma_PhaseI()
{
	// vvvv  build a feasible gamma, s.t. P*z_warm - h <= gamma holds, choose gamma = (P*z_warm-h)+difference
	Matrix<Type, _m, 1> c_tmp;		// use this to temporarily copy c-blocks out of z
	Matrix<Type, _n, 1> x_bar_tmp;	// use this to temporarily copy x_bar-blocks out of z
	Matrix<Type, Dynamic, 1> diff;		// length of vector not yet fixed
	
	int tmp = 0;		// variable used to count position
	// special treatment at beginning	// c_tmp = z.segment(0,_m);
	c_tmp = z.template segment<_m>(0);
	diff = (fu[0] - Fu[0]*K*x_hat) - (Fu[0]*c_tmp);	// should be >0
	// cout << "diff:" << endl << diff << endl << endl;
	for (int j=1; j <= _nInp; j++)
	{
		gamma[j-1] = -diff[j-1] + difference;
	}
	tmp = tmp + _nInp;
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		x_bar_tmp = z.template segment<_n>((i-1)*offset+_m+_n);
		c_tmp = z.template segment<_m>((i-1)*offset+_m+_n+_n);		
		diff = fx[i-1] - Fx[i-1]*x_bar_tmp;
		for (int j=1; j<= _nSt; j++)
		{
			gamma[tmp+j-1] = -diff[j-1] + difference;
		}
		tmp = tmp + _nSt;		
		diff = fu[i] - (Fu[i]*K*x_bar_tmp + Fu[i]*c_tmp);
		for (int j=1; j <= _nInp; j++)
		{
			gamma[tmp+j-1] = -diff[j-1] + difference;
		}
		tmp = tmp + _nInp;
	}
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	diff = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	for (int j=1; j<= _nSt; j++)
	{
		gamma[tmp+j-1] = -diff[j-1] + difference;
	}
	tmp = tmp + _nSt;
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);	// depends on position of invariant set
	diff = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		gamma[tmp+j-1] = -diff[j-1] + difference;
	}
	// cout << setprecision(15) << "gamma:" << endl << gamma << endl << endl;
	// ^^^ built a gamma
}

// ------------ function computes primal and dual residua needed by phaseI ---------------------
// stores in r_p (same as in PhaseII) and r_d_hat
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compRdRp_PhaseI()
{
		// r_p = C*z - b
		r_p.template segment<_n>(0) = -Bm*z_hat.template segment<_m>(0) + z_hat.template segment<_n>(_m+_nInp) - ( Am_tilde + Bm_bar)*x_hat - tm_tilde;
		r_p.template segment<_n>(_n) = -B*z_hat.template segment<_m>(0) + z_hat.template segment<_n>(_m+_nInp+_n) - (A_bar*x_hat + s);
		int tmp = _nInp + _nSt;
		for (int i=1; i <= _N-1; i++)
		{
			r_p.template segment<_n>(2*_n*i) = -Am_tilde*z_hat.template segment<_n>((i-1)*offset+_m+tmp-_nSt) - Bm_bar*z_hat.template segment<_n>((i-1)*offset+_m+_n+tmp-_nSt) -
				Bm*z_hat.template segment<_m>((i-1)*offset+_m+_n+_n+tmp) + z_hat.template segment<_n>((i-1)*offset+_m+_n+_n+_m+tmp+_nInp) - tm_tilde;

			r_p.template segment<_n>(2*_n*i+_n) = -A_bar*z_hat.template segment<_n>((i-1)*offset+_m+_n+tmp-_nSt) - B*z_hat.template segment<_m>((i-1)*offset+_m+_n+_n+tmp) + 
							z_hat.template segment<_n>((i-1)*offset+_m+_n+_n+_m+_n+tmp+_nInp) - s;
			tmp = tmp + _nInp + _nSt;
		}
		// cout << setprecision(15) << "z_hat" << endl << z_hat << endl << endl;
		// cout << setprecision(15) << "r_p in PhaseI:" << endl << r_p << endl << endl;	
		
		
		// r_d_hat = 2*H_hat*z_hat + g_hat + kappa*P_hat'*d + C_hat'*nu
		
		// note: H_hat = reg_hat * Identity();
		// r_d_hat.template segment<_m+_nInp>(0) = 2*reg_hat* Matrix<Type, Dynamic, Dynamic, 0, _m+_nInp, _m+_nInp>::Identity(_m+_nInp, _m+_nInp)*z_hat.template segment<_m+_nInp>(0) + g_hat_c[0] + 
		// cout << "reg_hat: " << reg_hat << endl;
		// cout << "kappa: " << kappa << endl;
		// cout << "z_hat: " << endl << z_hat << endl << endl;
		// cout << "nu: " << endl << nu << endl << endl;
		
		r_d_hat.template segment<_m+_nInp>(0) = 2*reg_hat*z_hat.template segment<_m+_nInp>(0) + g_hat_c[0] + 
				kappa * Fu_hat_transp[0] * d.template segment<_nInp>(0) - Bm_hat_transp[0] * nu.template segment<_n>(0) - B_hat_transp[0] * nu.template segment<_n>(_n);
		int offset1 = 2*_n;	// offset required in nu for C'*nu
		tmp = _nInp;	// offset required in d for P'*d	
		// cout << "bye 1" << endl;
		for (int i=1; i<= _N-1; i++)
		{
			// cout << "round " << i << endl << endl;
			// r_d_hat.template segment<_n>(_m+(i-1)*offset+tmp) = 2*reg_hat* Matrix<Type, _n, _n>::Identity(_n, _n) *z_hat.template segment<_n>(_m+(i-1)*offset+tmp) +
			r_d_hat.template segment<_n>(_m+(i-1)*offset+tmp) = 2*reg_hat*z_hat.template segment<_n>(_m+(i-1)*offset+tmp) +
					nu.template segment<_n>((i-1)*offset1) - Am_tilde_transp*nu.template segment<_n>((i-1)*offset1+_n+_n);
			
			if (i != _pos_omega)
			{
				// r_d_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = 2*reg_hat*Matrix<Type, _n+_nSt, _n+_nSt>::Identity(_n+_nSt, _n+_nSt)*z_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) + g_hat_x[i-1] +
				r_d_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = 2*reg_hat*z_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) + g_hat_x[i-1] +
						kappa*Fx_hat_transp[i-1] * d.template segment<_nSt>(tmp) + kappa*Fu_bar_hat_transp[i-1]*d.template segment<_nInp>(tmp+_nSt) + 
						Identity_hat_transp[i-1]*nu.template segment<_n>((i-1)*offset1+_n) - Bm_bar_hat_transp[i-1]*nu.template segment<_n>((i-1)*offset1+_n+_n) - A_bar_hat_transp[i-1]*nu.template segment<_n>((i-1)*offset1+_n+_n+_n);		
			}
			else	// i == _pos_omega
			{
				// r_d_hat.template segment<_n+_nSt>(_m+(_pos_omega-1)*offset+_n+tmp) = 2*reg_hat*Matrix<Type, _n+_nSt, _n+_nSt>::Identity(_n+_nSt, _n+_nSt)*z_hat.template segment<_n+_nSt>(_m+(_pos_omega-1)*offset+_n+tmp) + g_hat_x[_pos_omega-1] +
				r_d_hat.template segment<_n+_nSt>(_m+(_pos_omega-1)*offset+_n+tmp) = 2*reg_hat*z_hat.template segment<_n+_nSt>(_m+(_pos_omega-1)*offset+_n+tmp) + g_hat_x[_pos_omega-1] +
						kappa*Fx_hat_transp[_pos_omega-1] * d.template segment<_nSt>(tmp) + kappa*Fu_bar_hat_transp[_pos_omega-1]*d.template segment<_nInp>(tmp+_nSt) + 
						Identity_hat_transp[_pos_omega-1]*nu.template segment<_n>((_pos_omega-1)*offset1+_n) - Bm_bar_hat_transp[_pos_omega-1]*nu.template segment<_n>((_pos_omega-1)*offset1+_n+_n) - A_bar_hat_transp[_pos_omega-1]*nu.template segment<_n>((_pos_omega-1)*offset1+_n+_n+_n) 
						+ F_xTheta_hat_transp*d.template segment<_nF_xTheta>(_N*(_nSt + _nInp));
			}
			// cout << "r_d_hat in PhaseI:" << endl << r_d_hat << endl << endl;	
			tmp = tmp + _nSt;
			// r_d_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) = 2*reg_hat*Matrix<Type, Dynamic, Dynamic, 0, _m+_nInp, _m+_nInp>::Identity(_m+_nInp, _m+_nInp)*z_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) + g_hat_c[i] +
			r_d_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) = 2*reg_hat*z_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) + g_hat_c[i] +
					kappa*Fu_hat_transp[i]*d.template segment<_nInp>(tmp) - 
					Bm_hat_transp[i]*nu.template segment<_n>((i-1)*offset1+_n+_n) - B_hat_transp[i]*nu.template segment<_n>((i-1)*offset1+_n+_n+_n);		
			tmp = tmp + _nInp;
			// 				cout << "r_d_hat in PhaseI:" << endl << r_d_hat << endl << endl;	
			
		}
		// r_d_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = 2*reg_hat* Matrix<Type, _n, _n>::Identity(_n, _n)*z_hat.template segment<_n>(_m+(_N-1)*offset+tmp) + 
		r_d_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = 2*reg_hat*z_hat.template segment<_n>(_m+(_N-1)*offset+tmp) + 
					 			nu.template segment<_n>((_N-1)*offset1);
		if (_pos_omega == _N)
		{
			// r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = 2*reg_hat*Matrix<Type, _n+_nSt, _n+_nSt>::Identity(_n+_nSt, _n+_nSt)*z_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) + g_hat_x[_N-1] + 
			r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = 2*reg_hat*z_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) + g_hat_x[_N-1] + 
								kappa*Fx_hat_transp[_N-1]*d.template segment<_nSt>(tmp) + kappa*F_xTheta_hat_transp*d.template segment<_nF_xTheta>(tmp+_nSt) + 
								Identity_hat_transp[_N-1]*nu.template segment<_n>((_N-1)*offset1+_n);
		}
		else
		{
			// r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = 2*reg_hat*Matrix<Type, _n+_nSt, _n+_nSt>::Identity(_n+_nSt, _n+_nSt)*z_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) + g_hat_x[_N-1] + 
			r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = 2*reg_hat*z_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) + g_hat_x[_N-1] + 
								kappa*Fx_hat_transp[_N-1]*d.template segment<_nSt>(tmp) + 
								Identity_hat_transp[_N-1]*nu.template segment<_n>((_N-1)*offset1+_n);
		}
		
		tmp = tmp + _nSt;
		// r_d_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = 2*reg_hat*Matrix<Type, _m+_nF_xTheta, _m+_nF_xTheta>::Identity(_m+_nF_xTheta, _m+_nF_xTheta)*z_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp)  + g_hat_theta + 
		r_d_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = 2*reg_hat*z_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp)  + g_hat_theta + 
								kappa * F_theta_hat_transp*d.template segment<_nF_xTheta>(tmp);
		// cout << setprecision(15) << "r_d_hat:" << endl << r_d_hat << endl << endl;
}

// ---------------- computes the updated for z_hat and nu
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDz_hatDnu_hat_PhaseI()
{
	compPhi_hat_PhaseI();
	compPhi_hat_tilde_PhaseI();	// Phi_hat_tilde = Phi_hat^{-1}
	compY_hat_PhaseI();	// Y is stored in Y[][] from PhaseII
	compBeta_hat_PhaseI();	// beta_hat stored in beta from PhaseII
	// cout << setprecision(15) << "beta: " << endl << beta << endl << endl;
	compL();	// elements stores in L from PhaseII
	// ******* error *********
	compDnu();	// elements stored in dnu from PhaseII, Y*dnu = -beta_hat;
	// cout << setprecision(15) << "dnu in PhaseI" << endl << dnu << endl << endl;
	compDz_hat_PhaseI();
	// cout << setprecision(15) << "dz_hat in PhaseI" << endl << dz_hat << endl << endl;
}

// ------- computes Phi_hat = 2*H_hat + kappa*P'*diag(d)^2*P ------------------
// ------------ works ------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compPhi_hat_PhaseI()
{
	// ------------ first, partition diag(d) into (2N+1) diagonal matrices called d_diag[2N+1]
	// !!!!!!!!!!!!!!! what if _nF_xTheta is not big enough?
	
	Matrix<Type, _m+_nInp, _m+_nInp> eye_mnInp;
	eye_mnInp.setIdentity();
	Matrix<Type, _n+_nSt, _n+_nSt> eye_nnSt;
	eye_nnSt.setIdentity();
	Matrix<Type, _m+_nF_xTheta, _m+_nF_xTheta> eye_mTheta;
	eye_mTheta.setIdentity();
	
	DiagonalMatrix<Type, Dynamic> d_diag[2*_N+1];
	int tmp = 0;
	for (int i=0 ; i <= _N-1 ; i++)
	{
		d_diag[2*i].diagonal() = d.template segment<_nInp>(tmp);
		tmp = tmp + _nInp;
		d_diag[2*i+1].diagonal() = d.template segment<_nSt>(tmp);
		tmp = tmp + _nSt;
	}
	d_diag[2*_N].diagonal() = d.template segment<_nF_xTheta>(tmp);
	
	// ------- compute elements of Phi_hat, where Pi_tilde from PhaseII can be used
	// Omicron_hat[0] = 2*reg_hat*Matrix<Type, Dynamic, Dynamic, 0, _m+_nInp, _m+_nInp>::Identity(_m+_nInp, _m+_nInp) + kappa*Fu_hat_transp[0]*d_diag[0]*d_diag[0]*Fu_hat[0];
	Omicron_hat[0] = 2*reg_hat*eye_mnInp + kappa*Fu_hat_transp[0]*d_diag[0]*d_diag[0]*Fu_hat[0];
	// consider integrating d_diag directly into this loop
	for (int i=1; i <= _N-1; i++)
	{
		if (i != _pos_omega)
		{
			Rho_hat[i-1] = 2*reg_hat*eye_nnSt  + kappa * ( Fx_hat_transp[i-1]*d_diag[2*(i-1)+1]*d_diag[2*(i-1)+1]*Fx_hat[i-1]  + Fu_bar_hat_transp[i-1]*d_diag[2*i]*d_diag[2*i]*Fu_bar_hat[i-1] );
		}
		else	// i == _pos_omega
		{
			Rho_hat[_pos_omega-1] = 2*reg_hat*eye_nnSt  + kappa * ( Fx_hat_transp[_pos_omega-1]*d_diag[2*(_pos_omega-1)+1]*d_diag[2*(_pos_omega-1)+1]*Fx_hat[_pos_omega-1]  + Fu_bar_hat_transp[_pos_omega-1]*d_diag[2*_pos_omega]*d_diag[2*_pos_omega]*Fu_bar_hat[_pos_omega-1] 
				+ F_xTheta_hat_transp * d_diag[2*_N]*d_diag[2*_N]*F_xTheta_hat );
		}
		Sigma_hat[i-1] = kappa * ( Fu_bar_hat_transp[i-1]*d_diag[2*i]*d_diag[2*i]*Fu_hat[i] );
		Omicron_hat[i] = 2*reg_hat*eye_mnInp + kappa * ( Fu_hat_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu_hat[i] );
	}
	// special treatment for last block
	if (_pos_omega == _N)
	{
		Rho_hat[_N-1] = 2*reg_hat*eye_nnSt + kappa * ( Fx_hat_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx_hat[_N-1] + F_xTheta_hat_transp*d_diag[2*_N]*d_diag[2*_N]*F_xTheta_hat );
	}
	else
	{
		Rho_hat[_N-1] = 2*reg_hat*eye_nnSt + kappa * ( Fx_hat_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx_hat[_N-1] );
	}
	Sigma_hat[_N-1] = kappa * ( F_xTheta_hat_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta_hat );
	Omicron_hat[_N] = 2*reg_hat*eye_mTheta + kappa * ( F_theta_hat_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta_hat );

	/*
	for (int i = 0; i <= _N-1; i++)
	{
		cout << "Omicron_hat[" << i << "]" << endl << Omicron_hat[i] << endl << endl;
		cout << "Rho_hat[" << i << "]" << endl << Rho_hat[i] << endl << endl;
		cout << "Sigma_hat[" << i << "]" << endl << Sigma_hat[i] << endl << endl;
	}
	cout << "Omicron_hat[" << _N << "]" << endl << Omicron_hat[_N] << endl << endl;
	*/
}

template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compPhi_hat_tilde_PhaseI()
{
	// cout << "hallo 1" << endl << endl;
	LOmicron_diag[0].compute(Omicron_hat[0]);
	LOmicron_hat_diag_transp[0] = LOmicron_diag[0].matrixLLT().transpose();
	Matrix<Type, _n, _n> eye;
	eye.setIdentity();
	
	for (int i=1 ; i<= _N-1; i++)
	{
		// cout << "hello 2" << endl << endl;
		LPi_diag[i-1].compute(2*reg_hat*eye);
		LPi_hat_diag_transp[i-1] = LPi_diag[i-1].matrixLLT().transpose();
		LRho_diag[i-1].compute(Rho_hat[i-1]);
		LRho_hat_diag_transp[i-1] = LRho_diag[i-1].matrixLLT().transpose();		
		LSigma_hat_offDiag_transp[i-1] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(Sigma_hat[i-1]);		
		LSigma_hat_offDiag[i-1] = LSigma_hat_offDiag_transp[i-1].transpose();
		LOmicron_diag[i].compute(Omicron_hat[i]-LSigma_hat_offDiag[i-1]*LSigma_hat_offDiag_transp[i-1]);
		LOmicron_hat_diag_transp[i] = LOmicron_diag[i].matrixLLT().transpose();
		
		if (i == _pos_omega)
		{
			// cout << "entering if in for loop" << endl << endl;
			LLambda0_hat_transp = LRho_diag[_pos_omega-1].matrixLLT().triangularView<Lower>().solve(Sigma_hat[_N-1]);
			// cout << "Sigma_hat[_N-1].rows(): " << Sigma_hat[_N-1].rows() << endl << endl;
			// cout << "Sigma_hat[_N-1].cols(): " << Sigma_hat[_N-1].cols() << endl << endl;
			// cout << "LRho_diag[_pos_omega-1].rows(): " << LRho_diag[_pos_omega-1].rows() << endl << endl;
			// cout << "LRho_diag[_pos_omega-1].cols(): " << LRho_diag[_pos_omega-1].cols() << endl << endl;
			// cout << "hi 1" << endl << endl;
			LLambda0_hat = LLambda0_hat_transp.transpose();
			// cout << "LLamba0_hat.rows(): " << LLambda0_hat.rows() << endl << endl;
			// cout << "LLamba0_hat.cols(): " << LLambda0_hat.cols() << endl << endl;
			// cout << "hi 2" << endl << endl;
			// cout << "LSigma_hat_offDiag[_pos_omega-1].rows()" << LSigma_hat_offDiag[_pos_omega-1].rows() << endl << endl;
			// cout << "LSigma_hat_offDiag[_pos_omega-1].cols()" << LSigma_hat_offDiag[_pos_omega-1].cols() << endl << endl;
			// cout << "LLambda0_hat_transp.rows()" << LLambda0_hat_transp.rows() << endl << endl;
			// cout << "LLambda0_hat_transp.cols()" << LLambda0_hat_transp.cols() << endl << endl;
			LLambda1_hat_transp = LOmicron_diag[_pos_omega].matrixLLT().triangularView<Lower>().solve(Matrix<Type,_m+_nInp,_m+_nF_xTheta>::Zero() - LSigma_hat_offDiag[_pos_omega-1]*LLambda0_hat_transp);
			// cout << "hi 3" << endl << endl;
			LLambda1_hat = LLambda1_hat_transp.transpose();
			// cout << "exiting if in for loop" << endl << endl;
		}
	}
	
	// cout << "hello 3" << endl << endl;
	LPi_diag[_N-1].compute(2*reg_hat*eye);
	LPi_hat_diag_transp[_N-1] = LPi_diag[_N-1].matrixLLT().transpose();
	LRho_diag[_N-1].compute(Rho_hat[_N-1]);
	LRho_hat_diag_transp[_N-1] = LRho_diag[_N-1].matrixLLT().transpose();		
	
	// cout << "hello 4" << endl << endl;
	if(_N == _pos_omega)
	{
		LSigma_hat_offDiag_transp[_N-1] = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(Sigma_hat[_N-1]);		
		LSigma_hat_offDiag[_N-1] = LSigma_hat_offDiag_transp[_N-1].transpose();
		LOmicron_diag[_N].compute(Omicron_hat[_N]-LSigma_hat_offDiag[_N-1]*LSigma_hat_offDiag_transp[_N-1]);
		LOmicron_hat_diag_transp[_N] = LOmicron_diag[_N].matrixLLT().transpose();
	}
	else
	{
		LOmicron_diag[_N].compute(Omicron_hat[_N] - LLambda0_hat*LLambda0_hat_transp - LLambda1_hat*LLambda1_hat_transp  );
		LOmicron_hat_diag_transp[_N] = LOmicron_diag[_N].matrixLLT().transpose();
	}
	
	
	/*
	for (int i=0; i<= _N-2; i++)
	{
		cout << setprecision (5) << "LOmicron_diag[" << i << ']' << endl << LOmicron_diag[i].matrixLLT() << endl;
		cout << setprecision (5) << "LPi_diag[" << i << ']' << endl << LPi_diag[i].matrixLLT() << endl;
		cout << setprecision (5) << "LRho_diag[" << i << ']' << endl << LRho_diag[i].matrixLLT() << endl;
		cout << setprecision (5) << "LSigma_hat_offDiag[" << i << ']' << endl << LSigma_hat_offDiag[i] << endl;
	}
	cout << setprecision (5) << "LOmicron_diag[" << _N-1 << ']' << endl << LOmicron_diag[_N-1].matrixLLT() << endl;
	cout << setprecision (5) << "LPi_diag[" << _N-1 << ']' << endl << LPi_diag[_N-1].matrixLLT() << endl;
	cout << setprecision (5) << "LRho_diag[" << _N-1 << ']' << endl << LRho_diag[_N-1].matrixLLT() << endl;
	if(_N == _pos_omega)
	{
		cout << setprecision (5) << "LSigma_hat_offDiag[" << _N-1 << ']' << endl << LSigma_hat_offDiag[_N-1] << endl;
	}
	// cout << "hello Kitty" << endl << endl;
	cout << setprecision (5) << "LOmicron_diag[" << _N << ']' << endl << LOmicron_diag[_N].matrixLLT() << endl;
	// cout << "bye Kitty" << endl << endl;
	if (_pos_omega != _N)
	{
		// cout << "hello" << endl << endl;
		cout << "LLambda0_hat" << endl << LLambda0_hat << endl << endl;
		cout << "LLambda1_hat" << endl << LLambda1_hat << endl << endl;
	}
	*/
	
}

// ------- computes Y_hat = C_hat * Phi_hat_tilde * C_hat' ------------------------
// Y_hat uses Y[][] since the dimensions are the same as in PhaseII
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compY_hat_PhaseI()
{
	// computation of Y is done in three steps: Y_hat = C_hat * X_hat
	// 1. L_hat*U_hat = C_hat'
	// 2. L_hat'*X_hat = U_hat --> L_hat*L_hat'*X_hat = Phi_hat*X_hat = C_hat', i.e. X_hat = Phi_hat_tilde * C_hat'
	// 3. C_hat*X_hat = C_hat*Phi_hat_tilde*C_hat' = Y_hat
	
	Matrix<Type, _n, _n> eye;
	eye.setIdentity();
	Matrix<Type, _m+_nInp, _n> zero_hat_mnInp_n;	// size must be adjusted
	zero_hat_mnInp_n.setZero();
	Matrix<Type, _m+_nF_xTheta,_n> zero_hat_mTheta_n;
	zero_hat_mTheta_n.setZero();
	
	// 1. Compute elements of Matrix U_hat
	// treat the first 2 U_bar and first 3 U specially
	U_hat[0] = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(-Bm_hat_transp[0]);
	U_bar_hat[0] = LPi_diag[0].matrixLLT().triangularView<Lower>().solve(eye);
	
	U_hat[1] = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(-B_hat_transp[0]);
	U_bar_hat[1] = LRho_diag[0].matrixLLT().triangularView<Lower>().solve(Identity_hat_transp[0]);
	// zero_hat = Matrix<Type, Dynamic, _n>::Zero(_m+_nInp,_n);
	U_hat[2] = LOmicron_diag[1].matrixLLT().triangularView<Lower>().solve( zero_hat_mnInp_n - LSigma_hat_offDiag[0]*U_bar_hat[1] );
	
	// remaining U_bar and U have structures
	for (int i=1; i<= _N-2; i++)
	{
		U_bar_hat[2+(i-1)*5] = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(-Am_tilde_transp);
		U_bar_hat[3+(i-1)*5] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(-Bm_bar_hat_transp[i-1]);
		U_hat[3+(i-1)*3] = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(-Bm_hat_transp[i] - LSigma_hat_offDiag[i-1]*U_bar_hat[3+(i-1)*5] );
		U_bar_hat[4+(i-1)*5] = LPi_diag[i].matrixLLT().triangularView<Lower>().solve( eye );
		
		U_bar_hat[5+(i-1)*5] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(-A_bar_hat_transp[i-1]);
		U_hat[4+(i-1)*3] = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(-B_hat_transp[i] - LSigma_hat_offDiag[i-1]*U_bar_hat[5+(i-1)*5]);
		U_bar_hat[6+(i-1)*5] = LRho_diag[i].matrixLLT().triangularView<Lower>().solve( Identity_hat_transp[i] );
		// zero_hat = Matrix<Type, Dynamic, _n>::Zero(_m+_nInp,_n);
		U_hat[5+(i-1)*3] = LOmicron_diag[i+1].matrixLLT().triangularView<Lower>().solve( zero_hat_mnInp_n - LSigma_hat_offDiag[i]*U_bar_hat[6+(i-1)*5] );
	}
	
	U_bar_hat[2+(_N-2)*5] = LPi_diag[_N-1-1].matrixLLT().triangularView<Lower>().solve(-Am_tilde_transp);
	U_bar_hat[3+(_N-2)*5] = LRho_diag[_N-1-1].matrixLLT().triangularView<Lower>().solve(-Bm_bar_hat_transp[_N-2]);
	U_hat[3+(_N-2)*3] = LOmicron_diag[_N-1].matrixLLT().triangularView<Lower>().solve(-Bm_hat_transp[_N-1] - LSigma_hat_offDiag[_N-1-1]*U_bar_hat[3+(_N-1-1)*5] );
	U_bar_hat[4+(_N-2)*5] = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve( eye );
	
	U_bar_hat[5+(_N-2)*5] = LRho_diag[_N-1-1].matrixLLT().triangularView<Lower>().solve(-A_bar_hat_transp[_N-1-1]);
	U_hat[4+(_N-2)*3] = LOmicron_diag[_N-1].matrixLLT().triangularView<Lower>().solve(-B_hat_transp[_N-1] - LSigma_hat_offDiag[_N-1-1]*U_bar_hat[5+(_N-1-1)*5]);
	U_bar_hat[6+(_N-2)*5] = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve( Identity_hat_transp[_N-1] );
	
	if (_N == _pos_omega)
	{
		// zero_hat = Matrix<Type, Dynamic, _n>::Zero(_m+_nF_xTheta,_n);
		U_hat[5+(_N-2)*3] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( zero_hat_mTheta_n - LSigma_hat_offDiag[_N-1]*U_bar_hat[6+(_N-1-1)*5] );
	}
	else
	{
		// zero_hat = Matrix<Type, Dynamic, _n>::Zero(_m+_nF_xTheta,_n);
		UO_hat[0] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(zero_hat_mTheta_n - LLambda0_hat*U_bar_hat[1+(_pos_omega-1)*5] - LLambda1_hat*U_hat[2+(_pos_omega-1)*3]);
		UO_hat[1] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(zero_hat_mTheta_n - LLambda0_hat*U_bar_hat[3+(_pos_omega-1)*5] - LLambda1_hat*U_hat[3+(_pos_omega-1)*3]);
		UO_hat[2] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(zero_hat_mTheta_n - LLambda0_hat*U_bar_hat[5+(_pos_omega-1)*5] - LLambda1_hat*U_hat[4+(_pos_omega-1)*3]);
	}
	
	/*
	cout << "U_hat[0]" << endl << U_hat[0] << endl << endl;
	cout << "U_bar_hat[0]" << endl << U_bar_hat[0] << endl << endl;
	cout << "U_hat[1]" << endl << U_hat[1] << endl << endl;
	cout << "U_bar_hat[1]" << endl << U_bar_hat[1] << endl << endl;
	cout << "U_hat[2]" << endl << U_hat[2] << endl << endl;
	for (int i=1; i<= _N-2; i++)
	{
		cout << "U_bar_hat[2+(i-1)*5]" << endl << U_bar_hat[2+(i-1)*5] << endl << endl ;
		cout << "U_bar_hat[3+(i-1)*5]" << endl << U_bar_hat[3+(i-1)*5] << endl << endl ;
		cout << "U_hat[3+(i-1)*3]" << endl << U_hat[3+(i-1)*3] << endl << endl;
		cout << "U_bar_hat[4+(i-1)*5]" << endl << U_bar_hat[4+(i-1)*5] << endl << endl;
		
		cout << "U_bar_hat[5+(i-1)*5]" << endl << U_bar_hat[5+(i-1)*5] << endl << endl ;
		cout << "U_hat[4+(i-1)*3]" << endl << U_hat[4+(i-1)*3] << endl << endl ;
		cout << "U_bar_hat[6+(i-1)*5]" << endl << U_bar_hat[6+(i-1)*5] << endl << endl ;
		cout << "U_hat[5+(i-1)*3]" << endl << U_hat[5+(i-1)*3] << endl << endl ;
	}
	cout << "U_bar_hat[2+(_N-2)*5]" << endl << U_bar_hat[2+(_N-2)*5] << endl << endl ;
	cout << "U_bar_hat[3+(_N-2)*5]" << endl << U_bar_hat[3+(_N-2)*5] << endl << endl ;
	cout << "U_hat[3+(_N-2)*3]" << endl << U_hat[3+(_N-2)*3] << endl << endl;
	cout << "U_bar_hat[4+(_N-2)*5]" << endl << U_bar_hat[4+(_N-2)*5] << endl << endl;
	
	cout << "U_bar_hat[5+(_N-2)*5]" << endl << U_bar_hat[5+(_N-2)*5] << endl << endl ;
	cout << "U_hat[4+(_N-2)*3]" << endl << U_hat[4+(_N-2)*3] << endl << endl ;
	cout << "U_bar_hat[6+(_N-2)*5]" << endl << U_bar_hat[6+(_N-2)*5] << endl << endl ;
	
	if(_N == _pos_omega)
	{
		cout << "U_hat[5+(_N-2)*3]" << endl << U_hat[5+(_N-2)*3] << endl << endl ;
	}
	else
	{
		cout << "UO_hat[0]" << endl << UO_hat[0] << endl << endl;
		cout << "UO_hat[1]" << endl << UO_hat[1] << endl << endl;
		cout << "UO_hat[2]" << endl << UO_hat[2] << endl << endl;
	}
	*/


	// 2. Compute elements in Matrix X
	// treat the first 2 X_bar and first 3 X specially
	
	if(_pos_omega != _N)
	{
		XO_hat[0] = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(UO_hat[0]);
		XO_hat[1] = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(UO_hat[1]);
		XO_hat[2] = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(UO_hat[2]);
	}
	
	X_hat[0] = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(U_hat[0]);
	X_bar_hat[0] = LPi_hat_diag_transp[0].template triangularView<Upper>().solve(U_bar_hat[0]);
 	
	if (_pos_omega == 1)
	{
		X_hat[1] = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(U_hat[1]);
		X_hat[2] = LOmicron_hat_diag_transp[1].template triangularView<Upper>().solve(U_hat[2] - LLambda1_hat_transp*XO_hat[0]);
		X_bar_hat[1] = LRho_hat_diag_transp[0].template triangularView<Upper>().solve(U_bar_hat[1] - LSigma_hat_offDiag_transp[0]*X_hat[2] - LLambda0_hat_transp*XO_hat[0]);
	}
	else
	{
		X_hat[1] = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(U_hat[1]);
		X_hat[2] = LOmicron_hat_diag_transp[1].template triangularView<Upper>().solve(U_hat[2]);
		X_bar_hat[1] = LRho_hat_diag_transp[0].template triangularView<Upper>().solve(U_bar_hat[1] - LSigma_hat_offDiag_transp[0]*X_hat[2]);
	}
	
	// remaining X_bar and X have structures
	// if(_N != _pos_omega), then the off-diagonal element incluences three columns, called col1, col2, col3
	for (int i = 1; i <= _N-2; i++)
	{
		// cout << "round: " << i << endl;
		
		if (i == _pos_omega)	// col2
		{
			X_bar_hat[2+(_pos_omega-1)*5] = LPi_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar_hat[2+(_pos_omega-1)*5]);
			X_hat[3+(_pos_omega-1)*3] = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_hat[3+(_pos_omega-1)*3] - LLambda1_hat_transp*XO_hat[1]);
			X_bar_hat[3+(_pos_omega-1)*5] = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar_hat[3+(_pos_omega-1)*5] - LSigma_hat_offDiag_transp[_pos_omega-1]*X_hat[3+(_pos_omega-1)*3] - LLambda0_hat_transp*XO_hat[1]);
			X_bar_hat[4+(_pos_omega-1)*5] = LPi_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar_hat[4+(_pos_omega-1)*5]);
		}
		else // standard
		{
			X_bar_hat[2+(i-1)*5] = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(U_bar_hat[2+(i-1)*5]);
			X_hat[3+(i-1)*3] = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(U_hat[3+(i-1)*3]);
			X_bar_hat[3+(i-1)*5] = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(U_bar_hat[3+(i-1)*5] - LSigma_hat_offDiag_transp[i-1]*X_hat[3+(i-1)*3]);
			X_bar_hat[4+(i-1)*5] = LPi_hat_diag_transp[i].template triangularView<Upper>().solve(U_bar_hat[4+(i-1)*5]);
		}
		
		if (i == _pos_omega)	// col3
		{
			X_hat[4+(_pos_omega-1)*3] = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve( U_hat[4+(_pos_omega-1)*3] - LLambda1_hat_transp*XO_hat[2] );
			X_bar_hat[5+(_pos_omega-1)*5] = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar_hat[5+(_pos_omega-1)*5] - LSigma_hat_offDiag_transp[_pos_omega-1]*X_hat[4+(_pos_omega-1)*3] - LLambda0_hat_transp*XO_hat[2] );
			X_hat[5+(_pos_omega-1)*3] = LOmicron_hat_diag_transp[_pos_omega+1].template triangularView<Upper>().solve(U_hat[5+(_pos_omega-1)*3]);
			X_bar_hat[6+(_pos_omega-1)*5] = LRho_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_bar_hat[6+(_pos_omega-1)*5] - LSigma_hat_offDiag_transp[_pos_omega]*X_hat[5+(_pos_omega-1)*3]);
		}
		else if (i == _pos_omega - 1) // col1
		{
			X_hat[4+(_pos_omega-2)*3] = LOmicron_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve( U_hat[4+(_pos_omega-2)*3] );
			X_bar_hat[5+(_pos_omega-2)*5] = LRho_hat_diag_transp[_pos_omega-2].template triangularView<Upper>().solve(U_bar_hat[5+(_pos_omega-2)*5] - LSigma_hat_offDiag_transp[_pos_omega-2]*X_hat[4+(_pos_omega-2)*3]);
			X_hat[5+(_pos_omega-2)*3] = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(U_hat[5+(_pos_omega-2)*3] - LLambda1_hat_transp*XO_hat[0] );
			X_bar_hat[6+(_pos_omega-2)*5] = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(U_bar_hat[6+(_pos_omega-2)*5] - LSigma_hat_offDiag_transp[_pos_omega-1]*X_hat[5+(_pos_omega-2)*3] - LLambda0_hat_transp*XO_hat[0]);
		}
		else	// off-diag element of P*z<h has no influence
		{
			X_hat[4+(i-1)*3] = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve( U_hat[4+(i-1)*3] );
			X_bar_hat[5+(i-1)*5] = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(U_bar_hat[5+(i-1)*5] - LSigma_hat_offDiag_transp[i-1]*X_hat[4+(i-1)*3]);
			X_hat[5+(i-1)*3] = LOmicron_hat_diag_transp[i+1].template triangularView<Upper>().solve(U_hat[5+(i-1)*3]);
			X_bar_hat[6+(i-1)*5] = LRho_hat_diag_transp[i].template triangularView<Upper>().solve(U_bar_hat[6+(i-1)*5] - LSigma_hat_offDiag_transp[i]*X_hat[5+(i-1)*3]);
		}		
	}
	
	// compute last two columns
	if (_pos_omega == _N)
	{
		X_bar_hat[2+(_N-2)*5] = LPi_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[2+(_N-2)*5]);
		X_hat[3+(_N-2)*3] = LOmicron_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_hat[3+(_N-2)*3]);
		X_bar_hat[3+(_N-2)*5] = LRho_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[3+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-2]*X_hat[3+(_N-2)*3]);
		X_bar_hat[4+(_N-2)*5] = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar_hat[4+(_N-2)*5]);
		
		X_hat[4+(_N-2)*3] = LOmicron_hat_diag_transp[_N-1].template triangularView<Upper>().solve( U_hat[4+(_N-2)*3] );
		X_bar_hat[5+(_N-2)*5] = LRho_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[5+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-2]*X_hat[4+(_N-2)*3]);
		X_hat[5+(_N-2)*3] = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(U_hat[5+(_N-2)*3]);
		X_bar_hat[6+(_N-2)*5] = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar_hat[6+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-1]*X_hat[5+(_N-2)*3]);
	}
	else if (_pos_omega == _N-1)	// compute col2 and col3
	{
		X_bar_hat[2+(_N-2)*5] = LPi_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[2+(_N-2)*5]);
		X_hat[3+(_N-2)*3] = LOmicron_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_hat[3+(_N-2)*3] - LLambda1_hat_transp*XO_hat[1] );
		X_bar_hat[3+(_N-2)*5] = LRho_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[3+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-2]*X_hat[3+(_N-2)*3] - LLambda0_hat_transp*XO_hat[1] );
		X_bar_hat[4+(_N-2)*5] = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar_hat[4+(_N-2)*5]);
		
		X_hat[4+(_N-2)*3] = LOmicron_hat_diag_transp[_N-1].template triangularView<Upper>().solve( U_hat[4+(_N-2)*3] - LLambda1_hat_transp*XO_hat[2]);
		X_bar_hat[5+(_N-2)*5] = LRho_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[5+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-2]*X_hat[4+(_N-2)*3] - LLambda0_hat_transp*XO_hat[2]);
		X_bar_hat[6+(_N-2)*5] = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar_hat[6+(_N-2)*5]);
	}
	else 	// standard
	{
		X_bar_hat[2+(_N-2)*5] = LPi_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[2+(_N-2)*5]);
		X_hat[3+(_N-2)*3] = LOmicron_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_hat[3+(_N-2)*3] );
		X_bar_hat[3+(_N-2)*5] = LRho_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[3+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-2]*X_hat[3+(_N-2)*3]);
		X_bar_hat[4+(_N-2)*5] = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar_hat[4+(_N-2)*5]);
		
		X_hat[4+(_N-2)*3] = LOmicron_hat_diag_transp[_N-1].template triangularView<Upper>().solve( U_hat[4+(_N-2)*3]);
		X_bar_hat[5+(_N-2)*5] = LRho_hat_diag_transp[_N-2].template triangularView<Upper>().solve(U_bar_hat[5+(_N-2)*5] - LSigma_hat_offDiag_transp[_N-2]*X_hat[4+(_N-2)*3]);
		X_bar_hat[6+(_N-2)*5] = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(U_bar_hat[6+(_N-2)*5]);
	}
	
	/*
	cout << "X_hat[0]" << endl << X_hat[0] << endl << endl;
	cout << "X_bar_hat[0]" << endl << X_bar_hat[0] << endl << endl;
	cout << "X_hat[1]" << endl << X_hat[1] << endl << endl;
	cout << "X_bar_hat[1]" << endl << X_bar_hat[1] << endl << endl;
	cout << "X_hat[2]" << endl << X_hat[2] << endl << endl;
	for (int i=1; i<= _N-2; i++)
	{
		cout << "X_bar_hat[2+(i-1)*5]" << endl << X_bar_hat[2+(i-1)*5] << endl << endl ;
		cout << "X_bar_hat[3+(i-1)*5]" << endl << X_bar_hat[3+(i-1)*5] << endl << endl ;
		cout << "X_hat[3+(i-1)*3]" << endl << X_hat[3+(i-1)*3] << endl << endl;
		cout << "X_bar_hat[4+(i-1)*5]" << endl << X_bar_hat[4+(i-1)*5] << endl << endl;
		
		cout << "X_bar_hat[5+(i-1)*5]" << endl << X_bar_hat[5+(i-1)*5] << endl << endl ;
		cout << "X_hat[4+(i-1)*3]" << endl << X_hat[4+(i-1)*3] << endl << endl ;
		cout << "X_bar_hat[6+(i-1)*5]" << endl << X_bar_hat[6+(i-1)*5] << endl << endl ;
		cout << "X_hat[5+(i-1)*3]" << endl << X_hat[5+(i-1)*3] << endl << endl ;
	}
	cout << "X_bar_hat[2+(_N-2)*5]" << endl << X_bar_hat[2+(_N-2)*5] << endl << endl ;
	cout << "X_bar_hat[3+(_N-2)*5]" << endl << X_bar_hat[3+(_N-2)*5] << endl << endl ;
	cout << "X_hat[3+(_N-2)*3]" << endl << X_hat[3+(_N-2)*3] << endl << endl;
	cout << "X_bar_hat[4+(_N-2)*5]" << endl << X_bar_hat[4+(_N-2)*5] << endl << endl;
	
	cout << "X_bar_hat[5+(_N-2)*5]" << endl << X_bar_hat[5+(_N-2)*5] << endl << endl ;
	cout << "X_hat[4+(_N-2)*3]" << endl << X_hat[4+(_N-2)*3] << endl << endl ;
	cout << "X_bar_hat[6+(_N-2)*5]" << endl << X_bar_hat[6+(_N-2)*5] << endl << endl ;
	if (_pos_omega == _N)
	{
		cout << "X_hat[5+(_N-2)*3]" << endl << X_hat[5+(_N-2)*3] << endl << endl ;
	}
	else
	{
		cout << "XO_hat[0]" << endl << XO_hat[0] << endl << endl;
		cout << "XO_hat[1]" << endl << XO_hat[1] << endl << endl;
		cout << "XO_hat[2]" << endl << XO_hat[2] << endl << endl;
	}
	*/	
	
	// 3. Compute Y = C*X
	// compute first three Y separately
	Y[0][0] = -Bm_hat[0]*X_hat[0] + X_bar_hat[0];
	
	Y[0][1] = -Bm_hat[0]*X_hat[1];
	Y[1][1] = -B_hat[0]*X_hat[1] + Identity_hat[0]*X_bar_hat[1];
	// compute rest by filling Y column by column; done for 2 neighboring rows
	// we fill Y column by column, treating the first two columns specially
	for (int i=1; i <= _N-1; i++)
	{
		Y[0][2*(i-1)+2] = X_bar_hat[2+(i-1)*5];
		Y[1][2*(i-1)+2] = Identity_hat[i-1]*X_bar_hat[3+(i-1)*5];
		Y[2][2*(i-1)+2] = -Am_tilde*X_bar_hat[2+(i-1)*5] - Bm_bar_hat[i-1]*X_bar_hat[3+(i-1)*5] - Bm_hat[i]*X_hat[3+(i-1)*3] + X_bar_hat[4+(i-1)*5];
			
		Y[0][2*(i-1)+3] = Identity_hat[i-1]*X_bar_hat[5+(i-1)*5];
		Y[1][2*(i-1)+3] = -Bm_bar_hat[i-1]*X_bar_hat[5+(i-1)*5] - Bm_hat[i]*X_hat[4+(i-1)*3];
		Y[2][2*(i-1)+3] = -A_bar_hat[i-1]*X_bar_hat[5+(i-1)*5] - B_hat[i]*X_hat[4+(i-1)*3] + Identity_hat[i]*X_bar_hat[6+(i-1)*5];
	}
	
	/*
	cout << "Y[0][0]" << endl << Y[0][0] << endl << endl; 
	cout << "Y[0][1]" << endl << Y[0][1] << endl << endl;
	cout << "Y[1][1]" << endl << Y[1][1] << endl << endl;;
	
	
	for (int i=1; i <= _N-1; i++)
	{
		cout << "Y[0][2*(i-1)+2]" << endl << Y[0][2*(i-1)+2] << endl << endl; 
		cout << " Y[1][2*(i-1)+2]" << endl << Y[1][2*(i-1)+2] << endl << endl;
		cout << "Y[2][2*(i-1)+2]" << endl << Y[2][2*(i-1)+2] << endl << endl;
			
		cout << "Y[0][2*(i-1)+3]" << endl << Y[0][2*(i-1)+3] << endl << endl; 
		cout << "Y[1][2*(i-1)+3]" << endl << Y[1][2*(i-1)+3] << endl << endl; 
		cout << "Y[2][2*(i-1)+3]" << endl << Y[2][2*(i-1)+3] << endl << endl; 
	}
	*/
}


// ------------ function computes beta_hat = -r_p_hat + C_hat*Phi_hat_tilde*r_d_hat --------
// beta_hat uses beta as variable since is same size
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compBeta_hat_PhaseI()
{
	// solve in steps
	// beta_hat = -r_p_hat + C_hat*tmp1_hat, where tmp1_hat: Phi_hat*tmp1_hat = r_d_hat
	// 1. L_hat * tmp2_hat = r_d_hat   --> compute tmp2_hat
	// 2. L_hat' * tmp1_hat = tmp2_hat  ---> compute tmp1_hat
	// 3. beta = -r_p + C_hat * tmp1_hat
	
	// 1. compute tmp2_hat: L_hat * tmp2_hat = r_d_hat
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> tmp2_hat;	// long vector
	// tmp2_hat.resize((_N*(_m + _n + _n) + _m) + num_constr);	
	tmp2_hat.template segment<_m+_nInp>(0) = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+_nInp>(0) );
	int tmp = _nInp;
	for (int i = 1; i <= _N-1; i++)
	{
		tmp2_hat.template segment<_n>(_m+(i-1)*offset+tmp) = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n>(_m+(i-1)*offset+tmp));
		tmp2_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp));
		tmp = tmp + _nSt;
		tmp2_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) - LSigma_hat_offDiag[i-1]*tmp2_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp-_nSt));
		tmp = tmp + _nInp;
	}
	if(_pos_omega == _N)
	{	
		tmp2_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		tmp2_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp));
		tmp = tmp + _nSt;
		tmp2_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) - LSigma_hat_offDiag[_N-1]*tmp2_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp-_nSt));
	}
	else	// use LLambda0_hat and LLambda1_hat
	{
		// cout << "hello1" << endl << endl;
		tmp2_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		// cout << "hello2" << endl << endl;
		tmp2_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp));
		// cout << "hello3" << endl << endl;
		tmp = tmp + _nSt;
		// cout << "hello 4" << endl << endl;
		tmp2_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) - LLambda0_hat*tmp2_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n)    - LLambda1_hat*tmp2_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt)   );
		// cout << "hello 5" << endl << endl;
	}
	// cout << "tmp2_hat:" << endl << tmp2_hat << endl << endl;
	
	// 2. compute tmp1: L'*tmp1 = tmp2
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> tmp1_hat;	// long vector
	// tmp1_hat.resize((_N*(_m + _n + _n) + _m) + num_constr);
	tmp1_hat.template segment<_m+_nInp>(0) = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+_nInp>(0));
	tmp = _nInp;
	for (int i = 1; i <= _pos_omega-1; i++)
	{
		tmp1_hat.template segment<_n>(_m+(i-1)*offset+tmp) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+(i-1)*offset+tmp));
		tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) - LSigma_hat_offDiag_transp[i-1]*tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		tmp = tmp + _nSt + _nInp;
	}
	if (_pos_omega != _N)	// means that we skipped a part
	{
		tmp = tmp + _nSt + _nInp;
	}
	// the missing block is computed after the last block is computed
	for (int i = _pos_omega+1; i <= _N-1; i++)
	{
		tmp1_hat.template segment<_n>(_m+(i-1)*offset+tmp) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+(i-1)*offset+tmp));
		tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) - LSigma_hat_offDiag_transp[i-1]*tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		tmp = tmp + _nSt + _nInp;
	}
	if (_pos_omega == _N)
	{
		tmp1_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt));
		tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) - LSigma_hat_offDiag_transp[_N-1]*tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt));
	}
	else	// standard ending, compute the missing blocks
	{
		tmp1_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt));
		tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp));
		
		tmp1_hat.template segment<_n>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)) = LPi_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)) );
		tmp1_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) - LLambda1_hat_transp*tmp1_hat.template segment<_m+_nF_xTheta>(_m+_nInp+(_N-1)*(offset+_nSt+_nInp)+_n+_n+_nSt)  );
		tmp1_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n) = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n) - LSigma_hat_offDiag_transp[_pos_omega-1]*tmp1_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) - LLambda0_hat_transp*tmp1_hat.template segment<_m+_nF_xTheta>(_m+_nInp+(_N-1)*(offset+_nSt+_nInp)+_n+_n+_nSt)  );
	}
	// cout << "tmp1_hat:" << endl << tmp1_hat << endl << endl;
	
	
	// 3. beta = -r_p + C_hat * tmp1_hat
	beta.template segment<_n>(0) = -r_p.template segment<_n>(0) + ( -Bm_hat[0]*tmp1_hat.template segment<_m+_nInp>(0) + tmp1_hat.template segment<_n>(_m+_nInp) );
	beta.template segment<_n>(_n) = -r_p.template segment<_n>(_n) + (- B_hat[0]*tmp1_hat.template segment<_m+_nInp>(0) + Identity_hat[0]*tmp1_hat.template segment<_n+_nSt>(_m+_n+_nInp) );
	tmp = _nInp;
	for(int i=1; i<= _N-1; i++)
	{
		beta.template segment<_n>(2*i*_n) = -r_p.template segment<_n>(2*i*_n) + (    
			- Am_tilde*tmp1_hat.template segment<_n>(_m+(i-1)*offset+tmp) 
			- Bm_bar_hat[i-1] * tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp)
			- Bm_hat[i]* tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt) 
			+ tmp1_hat.template segment<_n>(i*offset+_m+tmp+_nInp+_nSt)    );
					
		beta.template segment<_n>(2*i*_n+_n) = -r_p.template segment<_n>( 2*i*_n+_n) + (  
		    - A_bar_hat[i-1] * tmp1_hat.template segment<_n+_nSt>(_m+_n+(i-1)*offset+tmp)
		    - B_hat[i] * tmp1_hat.template segment<_m+_nInp>(_m+_n+(i-1)*offset+_n+tmp+_nSt)
			+ Identity_hat[i]*tmp1_hat.template segment<_n+_nSt>(_m+(i)*offset+_n+tmp+_nInp+_nSt)        );
		tmp = tmp + _nSt + _nInp;
	}
	// cout << setprecision(15) << "beta PhaseI:" << endl << beta << endl << endl;
}


// ---- function computes Phi_hat * dz_hat = -r_d_hat - C_hat' * dnu ----------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDz_hat_PhaseI()
{
	// computed in parts
	// 1. tmp_hat = -r_d_hat - C_hat' * dnu
	// 2. L*L'*dz_hat = tmp_hat
	// 3. L*tmp1_hat = tmp_hat
	// 4. L'*dz_hat = tmp1_hat
	
	// 1. compute tmp_hat
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> tmp_hat;	// long vector
	// tmp_hat.resize(_N*(_m + _n + _n) + _m + num_constr);
	
	tmp_hat.template segment<_m+_nInp>(0) = -r_d_hat.template segment<_m+_nInp>(0) + Bm_hat_transp[0]*dnu.template segment<_n>(0) + B_hat_transp[0]*dnu.template segment<_n>(_n);
	int tmp = _nInp;
	for (int i=1; i<= _N-1; i++)
	{
		tmp_hat.template segment<_n>(_m+(i-1)*offset+tmp) = -r_d_hat.template segment<_n>(_m+(i-1)*offset+tmp) - dnu.template segment<_n>(2*(i-1)*_n) + Am_tilde_transp*dnu.template segment<_n>(2*(i-1)*_n+_n+_n) ;
		tmp_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = -r_d_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) - Identity_hat_transp[i-1]*dnu.template segment<_n>(2*(i-1)*_n+_n) + Bm_bar_hat_transp[i-1]*dnu.template segment<_n>(2*(i-1)*_n+_n+_n) + A_bar_hat_transp[i-1] * dnu.template segment<_n>(2*(i-1)*_n+_n+_n+_n);
		tmp = tmp + _nSt;
		tmp_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) = -r_d_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) + Bm_hat_transp[i]*dnu.template segment<_n>(2*(i-1)*_n+_n+_n) + B_hat_transp[i]*dnu.template segment<_n>(2*(i-1)*_n+_n+_n+_n);
		tmp = tmp + _nInp;
	}
	tmp_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = -r_d_hat.template segment<_n>(_m+(_N-1)*offset+tmp) - dnu.template segment<_n>(2*(_N-1)*_n);
	tmp_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = -r_d_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) - Identity_hat_transp[_N-1]*dnu.template segment<_n>(2*(_N-1)*_n+_n);
	tmp = tmp + _nSt;
	tmp_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = -r_d_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp);
	// cout << "tmp_hat:" << endl << tmp_hat << endl << endl;

	
	// 3. L*tmp1_hat = tmp_hat
	Matrix<Type, (_N*(_nInp+_nSt)+_nF_xTheta) + (_N*(_m + _n + _n) + _m) , 1> tmp1_hat;	// long vector
	// tmp1_hat.resize(_N*(_m + _n + _n) + _m + num_constr);
	tmp1_hat.template segment<_m+_nInp>(0) = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+_nInp>(0) );
	tmp = _nInp;
	for (int i = 1; i <= _N-1; i++)
	{
		tmp1_hat.template segment<_n>(_m+(i-1)*offset+tmp) = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n>(_m+(i-1)*offset+tmp));
		tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp));
		tmp = tmp + _nSt;
		tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp) - LSigma_hat_offDiag[i-1]*tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp-_nSt));
		tmp = tmp + _nInp;
	}
	if(_pos_omega == _N)
	{
		tmp1_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp));
		tmp = tmp + _nSt;
		tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) - LSigma_hat_offDiag[_N-1]*tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp-_nSt));
	}
	else
	{
		tmp1_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp));
		tmp = tmp + _nSt;
		tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp) - LLambda0_hat*tmp1_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n)  - LLambda1_hat*tmp1_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) );
	}
	
	// cout << "tmp1_hat:" << endl << tmp1_hat << endl << endl;
	
	
	// 4. L'*dz_hat = tmp1_hat
	dz_hat.template segment<_m+_nInp>(0) = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+_nInp>(0));
	tmp = _nInp;
	for (int i = 1; i <= _pos_omega-1; i++)
	{
		dz_hat.template segment<_n>(_m+(i-1)*offset+tmp) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+(i-1)*offset+tmp));
		dz_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		dz_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) - LSigma_hat_offDiag_transp[i-1]*dz_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		tmp = tmp + _nSt + _nInp;
	}
	if (_pos_omega != _N)	// means that we skipped a part
	{
		tmp = tmp + _nSt + _nInp;
	}
	for (int i = _pos_omega+1; i <= _N-1; i++)
	{
		dz_hat.template segment<_n>(_m+(i-1)*offset+tmp) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+(i-1)*offset+tmp));
		dz_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		dz_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+_nSt>(_m+(i-1)*offset+_n+tmp) - LSigma_hat_offDiag_transp[i-1]*dz_hat.template segment<_m+_nInp>(_m+(i-1)*offset+_n+_n+tmp+_nSt));
		tmp = tmp + _nSt + _nInp;
	}
	
	if (_pos_omega == _N)
	{
		dz_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		dz_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt));
		dz_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) - LSigma_hat_offDiag_transp[_N-1]*dz_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt));
	}
	else 	// standard ending, compute missing block in middle
	{
		dz_hat.template segment<_n>(_m+(_N-1)*offset+tmp) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+(_N-1)*offset+tmp));
		dz_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+_nF_xTheta>(_m+(_N-1)*offset+_n+_n+tmp+_nSt));
		dz_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+_nSt>(_m+(_N-1)*offset+_n+tmp));

		dz_hat.template segment<_n>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)) = LPi_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)) );
		dz_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) - LLambda1_hat_transp*dz_hat.template segment<_m+_nF_xTheta>(_m+_nInp+(_N-1)*(offset+_nSt+_nInp)+_n+_n+_nSt)  );
		dz_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n) = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n) - LSigma_hat_offDiag_transp[_pos_omega-1]*dz_hat.template segment<_m+_nInp>(_m+_nInp+(_pos_omega-1)*(offset+_nSt+_nInp)+_n+_n+_nSt) - LLambda0_hat_transp*dz_hat.template segment<_m+_nF_xTheta>(_m+_nInp+(_N-1)*(offset+_nSt+_nInp)+_n+_n+_nSt)  );
	}
	
	// cout << "dz_hat:" << endl << dz_hat << endl << endl;
}

// ------------ function computes a feasibe gamma needed by phaseI ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: testPos_PhaseI()
{
	Matrix<Type, Dynamic, 1> c_hat_tmp;		// use this to temporarily copy c-blocks out of z
	Matrix<Type, Dynamic, 1, 0, _n+_nSt, 1> x_bar_hat_tmp;	// use this to temporarily copy x_bar-blocks out of z
	Matrix<Type, Dynamic, 1> check;		// length of vector not yet fixed
	int tmp = 0;		// variable used to count position and offset
	
	// special treatment at beginning
	c_hat_tmp = z_hat.template segment<_m+_nInp>(0);
	check = (fu[0] - Fu[0]*K*x_hat) - (Fu_hat[0]*c_hat_tmp);	// should be >0
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j <= _nInp; j++)
	{
		if (check[j-1] <= 0)
		{
			return 0;
		}
	}
	tmp = tmp + _nInp;
	
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		// cout << "entering round: " << i << " in testPos_PhaseI()" << endl << endl;
		x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((i-1)*offset+_m+_n+tmp);
		c_hat_tmp = z_hat.template segment<_m+_nInp>((i-1)*offset+_m+_n+_n+tmp+_nSt);		
		check = fx[i-1] - Fx_hat[i-1]*x_bar_hat_tmp;
		// cout << "check:" << endl << check << endl << endl;
		for (int j=1; j<= _nSt; j++)
		{
			if (check[j-1] <= 0)
			{
				return 0;
			}
		}
		check = fu[i] - (Fu_bar_hat[i-1]*x_bar_hat_tmp + Fu_hat[i]*c_hat_tmp);
		// cout << "check:" << endl << check << endl << endl;
		for (int j=1; j<=_nInp; j++)
		{
			if (check[j-1] <= 0)
			{
				return 0;
			}
		}
		tmp = tmp + _nSt + _nInp;	
	}
	// cout << "finished general treatment at the beginning" << endl << endl;
	
	x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((_N-1)*offset+_m+_n+tmp);
	c_hat_tmp = z_hat.template segment<_m+_nF_xTheta>((_N-1)*offset+_m+_n+_n+tmp+_nSt);
	check = fx[_N-1] - (Fx_hat[_N-1]*x_bar_hat_tmp);
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j<= _nSt; j++)
	{
		if (check[j-1] <= 0)
		{
			return 0;
		}
	}
	x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((_pos_omega-1)*(offset+_nInp+_nSt)+_m+_nInp+_n);
	check = f_xTheta - (F_xTheta_hat*x_bar_hat_tmp + F_theta_hat*c_hat_tmp);
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if (check[j-1] <= 0)
		{
			return 0;
		}
	}
	return 1;
}

// ------------ function computes (d_hat)_i = (h - P_hat*z_hat)_i  --------
// elements are stored in d
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compD_hat_PhaseI()
{
	Matrix<Type, _m+_nInp, 1> c_hat_tmp_mnInp;		// use this to temporarily copy c-blocks out of z
	Matrix<Type, _m+_nF_xTheta, 1> c_hat_tmp_mTheta;
	Matrix<Type, _n+_nSt, 1> x_bar_hat_tmp;	// use this to temporarily copy x_bar-blocks out of z
	Matrix<Type, Dynamic, 1> diff;		// length of vector not yet fixed
	int tmp = 0;		// variable used to count position and offset
	
	// special treatment at beginning
	c_hat_tmp_mnInp = z_hat.template segment<_m+_nInp>(0);
	diff = (fu[0] - Fu[0]*K*x_hat) - (Fu_hat[0]*c_hat_tmp_mnInp);
	// cout << "diff:" << endl << diff << endl << endl;
	for (int j=1; j <= _nInp; j++)
	{
		d[j-1] = 1/diff[j-1];
	}
	tmp = tmp + _nInp;
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((i-1)*offset+_m+_n+tmp);
		c_hat_tmp_mnInp = z_hat.template segment<_m+_nInp>((i-1)*offset+_m+_n+_n+tmp+_nSt);		
		diff = fx[i-1] - Fx_hat[i-1]*x_bar_hat_tmp;
		for (int j=1; j<= _nSt; j++)
		{
			d[tmp+j-1] = 1/diff[j-1];
		}
		tmp = tmp  + _nSt;
		
		diff = fu[i] - (Fu_bar_hat[i-1]*x_bar_hat_tmp + Fu_hat[i]*c_hat_tmp_mnInp);
		for (int j=1; j<=_nInp; j++)
		{
			d[tmp+j-1] = 1/diff[j-1];
		}
		tmp = tmp + _nInp;	
	}
	
	x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((_N-1)*offset+_m+_n+tmp);
	c_hat_tmp_mTheta = z_hat.template segment<_m+_nF_xTheta>((_N-1)*offset+_m+_n+_n+tmp+_nSt);
	diff = fx[_N-1] - (Fx_hat[_N-1]*x_bar_hat_tmp);
	for (int j=1; j<= _nSt; j++)
	{
		d[tmp+j-1] = 1/diff[j-1];
	}
	tmp = tmp + _nSt;
	x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((_pos_omega-1)*(offset+_nInp+_nSt)+_m+_nInp+_n);
	diff = f_xTheta - (F_xTheta_hat*x_bar_hat_tmp + F_theta_hat*c_hat_tmp_mTheta);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		d[tmp+j-1] = 1/diff[j-1];
	}
}

// ------------ function checks if gamma < 0 and extracts z_warm  --------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: testNeg_PhaseI()
{
	Matrix<Type, Dynamic, 1, 0, _nSt + _nInp, 1> check;		// stores gamma segments
	int tmp = 0;		// variable used to count position and offset
	// general treatment in the middle; offset = _n + _n + _m;
	
	cout << "starting new testNeg_PhaseI()" << endl;
	
	for (int i=1; i <= _N; i++) // blocks in the middle
	{	
		check.resize(_nInp + _nSt);
		check << z_hat.template segment<_nInp>((i-1)*offset+_m+tmp) ,
		         z_hat.template segment<_nSt>((i-1)*offset+_m+tmp+_n+_n+_nInp);
		cout << "check in testNeg_PhaseI()" << endl << check << endl << endl;
		for (int j=1; j<= _nInp+_nSt; j++)
		{
			if (check[j-1] >= 0)
			{
				// cout << "z_warm: " << endl << z_warm << endl << endl;
				return 0;
			}
		}
		z_warm.template segment<_m>((i-1)*offset) = z_hat.template segment<_m>((i-1)*offset+tmp);
		z_warm.template segment<_n+_n>((i-1)*offset + _m) = z_hat.template segment<_n+_n>((i-1)*offset+_m+tmp+_nInp);
		tmp = tmp + _nSt + _nInp;	
	}
	
	check.resize(_nF_xTheta);
	check = z_hat.template segment<_nF_xTheta>(_N*offset + _m + tmp);
	cout << "check in the end testNeg_PhaseI()" << endl << check << endl << endl;
	// cout << "check in testNeg_PhaseI" << endl << check << endl << endl;
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if (check[j-1] >= 0)
		{
			// cout << "z_warm: " << endl << z_warm << endl << endl;
			return 0;
		}
	}
	z_warm.template segment<_m>(_N*offset ) = z_hat.template segment<_m>(_N*offset + tmp);
	return 1;
}


// ------------ function returns largest t: P_hat*(z_hat+t*dz_hat) < h ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
double LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compT_PhaseI()
{
	
	Matrix<Type, Dynamic, 1> c_hat_tmp;		// use this to temporarily copy c-blocks out of z_hat
	Matrix<Type, Dynamic, 1> dc_hat_tmp;	// copy blocks from dz_hat
	Matrix<Type, Dynamic, 1, 0, _n+_nSt, 1> x_bar_hat_tmp;	// use this to temporarily copy x_bar-blocks out of z_hat
	Matrix<Type, Dynamic, 1, 0, _n+_nSt, 1> dx_bar_hat_tmp;	// use this to temporarily copy x_bar-blocks out of dz_hat
	
	Matrix<Type, Dynamic, 1> check;		// check = h - P_hat*z_hat
	Matrix<Type, Dynamic, 1> dcheck;	// dcheck = P_hat*dz_hat
	Matrix<Type, Dynamic, 1> t_vec;		// t_vec = check ./ dcheck
	double t=1;		// stores the smallest variable
	double t_tmp;	// t_tmp = min(t_vec)
	
	int tmp = 0;		// variable used to count position and offset
	c_hat_tmp = z_hat.template segment<_m+_nInp>(0);
	dc_hat_tmp = dz_hat.template segment<_m+_nInp>(0);
	check = (fu[0] - Fu[0]*K*x_hat) - (Fu_hat[0]*c_hat_tmp);	// should be >0
	dcheck = Fu_hat[0]*dc_hat_tmp;
	t_vec.resize(_nInp);
	t_vec.setConstant(1);
	for (int j=1; j <= _nInp; j++)
	{
		if (dcheck[j-1] > 0)	// neg. cases not interesting
		{
			t_vec[j-1] = check[j-1]/dcheck[j-1];
		}
	}
	t_tmp = t_vec.minCoeff();
	if (t_tmp < t)
	{
		t = t_tmp;
	}
	tmp = tmp + _nInp;
	
	
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		// cout << "entering round: " << i << " in testPos_PhaseI()" << endl << endl;
		x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((i-1)*offset+_m+_n+tmp);
		dx_bar_hat_tmp = dz_hat.template segment<_n+_nSt>((i-1)*offset+_m+_n+tmp);
		c_hat_tmp = z_hat.template segment<_m+_nInp>((i-1)*offset+_m+_n+_n+tmp+_nSt);
		dc_hat_tmp = dz_hat.template segment<_m+_nInp>((i-1)*offset+_m+_n+_n+tmp+_nSt);		
		check = fx[i-1] - Fx_hat[i-1]*x_bar_hat_tmp;
		dcheck = Fx_hat[i-1]*dx_bar_hat_tmp;
		t_vec.resize(_nSt);
		t_vec.setConstant(1);
		// cout << "check:" << endl << check << endl << endl;
		for (int j=1; j<= _nSt; j++)
		{
			if(dcheck[j-1]>0)
			{
				t_vec[j-1] = check[j-1]/dcheck[j-1];
			}
		}
		t_tmp = t_vec.minCoeff();
		if (t_tmp < t)
		{
			t = t_tmp;
		}
		
		check = fu[i] - (Fu_bar_hat[i-1]*x_bar_hat_tmp + Fu_hat[i]*c_hat_tmp);
		dcheck = Fu_bar_hat[i-1]*dx_bar_hat_tmp + Fu_hat[i]*dc_hat_tmp;
		t_vec.resize(_nInp);
		t_vec.setConstant(1);
		for (int j=1; j<=_nInp; j++)
		{
			if(dcheck[j-1]>0)
			{
				t_vec[j-1] = check[j-1]/dcheck[j-1];
			}
		}
		t_tmp = t_vec.minCoeff();
		if (t_tmp < t)
		{
			t = t_tmp;
		}
		tmp = tmp + _nSt + _nInp;	
	}
	// cout << "finished general treatment at the beginning" << endl << endl;
	
	// special case for last blocks
	x_bar_hat_tmp = z_hat.template segment<_n+_nSt>((_N-1)*offset+_m+_n+tmp);
	dx_bar_hat_tmp = dz_hat.template segment<_n+_nSt>((_N-1)*offset+_m+_n+tmp);
	c_hat_tmp = z_hat.template segment<_m+_nF_xTheta>((_N-1)*offset+_m+_n+_n+tmp+_nSt);
	dc_hat_tmp = dz_hat.template segment<_m+_nF_xTheta>((_N-1)*offset+_m+_n+_n+tmp+_nSt);
	check = fx[_N-1] - (Fx_hat[_N-1]*x_bar_hat_tmp);
	dcheck = Fx_hat[_N-1]*dx_bar_hat_tmp;
	t_vec.resize(_nSt);
	t_vec.setConstant(1);
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j<= _nSt; j++)
	{
		if (dcheck[j-1]>0)
		{
			t_vec[j-1] = check[j-1]/dcheck[j-1];
		}
	}
	t_tmp = t_vec.minCoeff();
	if (t_tmp < t)
	{
		t = t_tmp;
	}
	
	x_bar_hat_tmp = z_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nInp+_nSt)+_n);
	dx_bar_hat_tmp = dz_hat.template segment<_n+_nSt>(_m+_nInp+(_pos_omega-1)*(offset+_nInp+_nSt)+_n);
	check = f_xTheta - (F_xTheta_hat*x_bar_hat_tmp + F_theta_hat*c_hat_tmp);
	dcheck = F_xTheta_hat*dx_bar_hat_tmp + F_theta_hat*dc_hat_tmp;
	
	t_vec.resize(_nF_xTheta);
	t_vec.setConstant(1);
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if(dcheck[j-1]>0)
		{
			t_vec[j-1] = check[j-1]/dcheck[j-1];
		}
	}
	t_tmp = t_vec.minCoeff();
	if (t_tmp < t)
	{
		t = t_tmp;
	}
	
	// return result
	if (t == 1)
	{
		return 1;
	} 
	else
	{
		return 0.99*t;	// guarantees strict feasibility
	}
}







#endif