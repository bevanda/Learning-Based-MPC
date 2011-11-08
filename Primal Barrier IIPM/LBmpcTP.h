// Ompctp.h -- a Learning Based MPC class template
// header for class LBmpcTP
// date: October 28, 2011
// author: Xiaojing ZHANG
// version: 0.1


/* Remarks:

2) horizon >= 3 assumed
3) issues with LLT, L_diag[2*Horizon] for both Y and Phi --> save it into normal matrices
4) issues of defining length and Type for vector d, use diagonal Matrix for d^2
5) Optimize index calculations, merge different loops into 1 loop
6) avoid using namespace Eigen/std
9) upper bound size of d_diag in compPhi and compPhi_hat_PhaseI()
10) omit isFeasible in step() and phaseI() and testPos() and testPos_hat();
11) compGamma(): use vector addition rather than piece by piece addition?
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
	int n_iter;		// number of Newton iterations for fixed kappa in PhaseII
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
	
	// !!!!!!!! Variables to represent Phi = L_Phi*L_Phi'
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
	
	Matrix<Type, _m, 1> dc_tmp;		// use this to temporarily copy c-blocks out of dz
	Matrix<Type, _n, 1> dx_bar_tmp;	// use this to temporarily copy x_bar-blocks out of dz
	Matrix<Type, _nSt, 1> dcheckX;	// dcheck = P*dz
	Matrix<Type, _nInp, 1> dcheckU;
	Matrix<Type, _nF_xTheta, 1> dcheckTheta;
	Matrix<Type, _nSt, 1> t_vecX;		// t_vec = check ./ dcheck
	Matrix<Type, _nInp, 1> t_vecU;
	Matrix<Type, _nF_xTheta, 1> t_vecTheta;
	Matrix<Type, _m+1, 1> dc_hat_tmp;	// copy blocks from dz_hat
	Matrix<Type, _n+1, 1> dx_bar_hat_tmp;	// use this to temporarily copy x_bar-blocks out of dz_hat
	
	
// ---------- vvv  Matrices for PhaseI() -------------
	double kappa_start_PhaseI;
	double weight_cost_PhaseI;
	int offset2;
	double weight_PhaseI;	// weight 
	int n_iter_PhaseI;	// max. number of Newton iterations for fixed kappa in PhaseI
	double reg_hat;		// regularisation term in Matrix H_hat
	
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> z_hat;	// long vector
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> z_hat_orig;	// long vector
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> dz_hat;	// long vector
	Matrix<Type, 2*_N*(_n+1)-1,1> nu_hat;
	Matrix<Type, 2*_N*(_n+1)-1,1> nu_hat_orig;
	Matrix<Type, 2*_N*(_n+1)-1,1> dnu_hat;

	Matrix<Type, 2*_N+1, 1> gamma;		// contains the s0, s1 , ..., s2N
	Matrix<Type, 2*_N*(_n+1)-1, 1> r_p_hat;
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> r_d_hat;	// long vector
	Matrix<Type, 2*_N*(_n+1)-1, 1> beta_hat;
	
	Matrix<Type, _nSt, _n+1> Fx_hat[_N];	// array of full-rank state constraint matrices
	Matrix<Type, _n+1, _nSt> Fx_hat_transp[_N];	// transpose of above
	Matrix<Type, _nInp, _m+1> Fu_hat[_N];	// array of (full-rank) input constraint matrices
	Matrix<Type, _m+1, _nInp> Fu_hat_transp[_N];	//transpose of above
	Matrix<Type, _nF_xTheta, _n+1> F_xTheta_hat;	// F_xTheta*x[m+N]+F_theta*theta <= f_xTheta
	Matrix<Type, _n+1, _nF_xTheta> F_xTheta_hat_transp;	// transpose of above
	Matrix<Type, _nF_xTheta, _m+1> F_theta_hat;	// full-rank constraint matrix on theta
	Matrix<Type, _m+1, _nF_xTheta> F_theta_hat_transp;	// transpose of above
	Matrix<Type, _nInp, _n+1> Fu_bar_hat[_N-1];	// array of (full-rank) input constraint matrices, #row may change, but <= nInp
	Matrix<Type, _n+1, _nInp> Fu_bar_hat_transp[_N-1];	//transpose of above
	
	Matrix<Type, _n+1, _m+1> Bm_hat;	// Bm_hat is not used for first block row in C_hat b/c is has a row less
	Matrix<Type, _m+1, _n+1> Bm_hat_transp;	// see above
	Matrix<Type, _n, _m+1> Bm_hat0;		// replaces the Bm_hat above for first block row
	Matrix<Type, _m+1, _n> Bm_hat0_transp;
	Matrix<Type, _n+1, _m+1> B_hat;
	Matrix<Type, _m+1, _n+1> B_hat_transp;
	Matrix<Type, _n+1, _n+1> Bm_bar_hat;
	Matrix<Type, _n+1, _n+1> Bm_bar_hat_transp;
	Matrix<Type, _n+1, _n+1> A_bar_hat;
	Matrix<Type, _n+1, _n+1> A_bar_hat_transp;
	Matrix<Type, _n+1, _n> Identity1_hat;		// Identity1_hat[0] is not used for first block row in C_hat b/c it has a row less
	Matrix<Type, _n, _n+1> Identity1_hat_transp;	// see above
	Matrix<Type, _n+1, _n+1> Identity2_hat;		
	Matrix<Type, _n+1, _n+1> Identity2_hat_transp;
	Matrix<Type, _n+1, _n> Am_tilde_hat;
	Matrix<Type, _n, _n+1> Am_tilde_hat_transp;
	
	Matrix<Type, _n+1, 1> x_bar_hat_tmp;	// temporary variables needed in compD_hat_PhaseI, compT_PhaseI, ...
	Matrix<Type, _m+1, 1> c_hat_tmp;		
	
	Matrix<Type, _m+1, _m+1> Omicron_hat[_N+1];
	Matrix<Type, _n+1, _n+1> Rho_hat[_N];
	Matrix<Type, _n+1, _m+1> Sigma_hat[_N];
	
	Matrix<Type, _m+1, _m+1> LOmicron_hat_diag_transp[_N+1];
	Matrix<Type, _n, _n> LPi_hat_diag_transp[_N];	// stores the transpose
	Matrix<Type, _n+1, _n+1> LRho_hat_diag_transp[_N];
	Matrix<Type, _m+1, _n+1> LSigma_hat_offDiag[_N];
	Matrix<Type, _n+1, _m+1> LSigma_hat_offDiag_transp[_N];
	Matrix<Type, _m+1, _n+1> LLambda0_hat;		// used for LLT decomposition if _pos_omega != _N
	Matrix<Type, _n+1, _m+1> LLambda0_hat_transp;
	Matrix<Type, _m+1, _m+1> LLambda1_hat;		// used for LLT decomposition if _pos_omega != _N
	Matrix<Type, _m+1, _m+1> LLambda1_hat_transp;
	
	Matrix<Type, Dynamic, Dynamic, 0, _n+1, _n+1> Y_hat[3][2*_N];	// Y[i][j] = Y_{i+1,j+1}, i=0,1,2
	Matrix<Type, Dynamic, Dynamic, 0, _n+1, _n+1> L_hat_diag_transp[2*_N];
	Matrix<Type, Dynamic, Dynamic, 0, _n+1, _n+1> L_hat_offDiag[2][2*_N];	// off-diag matrices are general square matrices
	Matrix<Type, Dynamic, Dynamic, 0, _n+1, _n+1> L_hat_offDiag_transp[2][2*_N];	// transpose of above
	
	Matrix<Type, Dynamic, Dynamic, 0, _m+1, _n+1> U_hat[3+(_N-1)*3];		// L*U = C', needed to compute Y, 
	Matrix<Type, Dynamic, Dynamic, 0, _n+1, _n+1> U_bar_hat[2+(_N-1)*5];	// L*U = C'
	Matrix<Type, _m+1, _n+1> UO_hat[3];
	Matrix<Type, Dynamic, Dynamic, 0, _m+1, _n+1> X_hat[3+(_N-1)*3];		// L'*X = U
	Matrix<Type, Dynamic, Dynamic, 0, _n+1, _n+1> X_bar_hat[2+(_N-1)*5];	// L'*X = U
	Matrix<Type, _m+1, _n+1> XO_hat[3];
	
	double difference;	// parameter needed to compute gamma
// ------------ ^^^ end Matrices PhaseI() --------------	
	
	


// ---------- private methods ----------
	// "testPos" tests if the given z_warm satisfies P*z_warm < h								
	bool testPos();		// returns 1, inequality is satisfied
						// return 0, if inequality is violated
	bool phaseI();	// in case P*z_warm >= h, phaseI will retrieve a good starting point
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
	void compL_hat();
	void compDnu_hat();
	void compDz_hat_PhaseI();
	// bool testPos_PhaseI();
	void compD_hat_PhaseI();
	bool testNeg_PhaseI();
	double compT_PhaseI();		// computes largest t: P_hat*(z_hat+t*dz_hat) < h
// --------- ^^^ functions for PhaseI() -------


	
public:
	LBmpcTP();	// default constructor
	
	// standard constructor: uses references whenever possible to avoid copying of big matrices
	// note: arrays cannot be copied by value, always given as pointers, length implicitly given by _N
	LBmpcTP(double kappa_arg, double kappa_PhaseI_arg, int n_iter_arg, int n_iter_PhaseI_arg, double mu_arg, double eps_barrier_arg, double eps_nt_arg, double eps_normRp_arg,
		  	double eps_ls_arg, double alpha_ls_arg, double beta_ls_arg, double reg_arg, double reg_PhaseI_arg, double weight_PhaseI_arg,
		   const Matrix<Type, _n, _n> &A_arg, Matrix<Type, _n, _m> &B_arg,
		   const Matrix<Type, _n, _n> &Q_tilde_arg, const Matrix<Type, _n, _n> &Q_tilde_f_arg,
		   const Matrix<Type, _m, _m> &R_arg, const Matrix<Type, _nSt, _n> Fx_arg[], 
		   const Matrix<Type, _nSt, 1> fx_arg[], const Matrix<Type, _nInp, _m> Fu_arg[],
		   const Matrix<Type, _nInp, 1> fu_arg[], const Matrix<Type, _nF_xTheta, _n> &F_xTheta_arg,
		   const Matrix<Type, _nF_xTheta, _m> &F_theta_arg, const Matrix<Type, _nF_xTheta, 1> &f_xTheta_arg,
		   const Matrix<Type, _m, _n> &K_arg, const Matrix<Type, _n, 1> &s_arg
	      );	
	

	// "step" computes and returns the optimal input 
	Matrix<Type, _m, 1> step(const Matrix<Type, _n, _n> &Lm_arg, const Matrix<Type, _n, _m> &Mm_arg, const Matrix<Type, _n, 1> &tm_arg,
														const Matrix<Type, _n, 1> &x_hat_arg,
														const Matrix<Type, _n, 1> x_star_arg[]
													   );
													
	//~LBmpcTP();	// destructor
};



//  ==================== Implementation of Methods ==================

// default constructor
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::LBmpcTP()
{
	// do nothing
}

// constructor
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::LBmpcTP(double kappa_arg, double kappa_PhaseI_arg, int n_iter_arg, int n_iter_PhaseI_arg,
	   double mu_arg, double eps_barrier_arg, double eps_nt_arg, double eps_normRp_arg, double eps_ls_arg, double alpha_ls_arg, double beta_ls_arg,
	   double reg_arg, double reg_PhaseI_arg, double weight_PhaseI_arg, const Matrix<Type, _n, _n> &A_arg, Matrix<Type, _n, _m> &B_arg,
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
	offset2 = offset+2;
	n_iter_PhaseI = n_iter_PhaseI_arg;
	weight_PhaseI = weight_PhaseI_arg;
	
	A = A_arg;
	A_transp = A_arg.transpose();
	B = B_arg;
	B_transp = B_arg.transpose();
	Q_tilde = Q_tilde_arg;
	Q_tilde_f = Q_tilde_f_arg;
	R = R_arg;	
	K = K_arg;
	K_transp = K.transpose();
		
	// num_constr = 0;
	for (int i=0; i < _N; i++)
	{
		Fx[i] = Fx_arg[i];
		Fx_transp[i] = Fx_arg[i].transpose();
		fx[i] = fx_arg[i];
		
		Fu[i] = Fu_arg[i];
		Fu_transp[i] = Fu_arg[i].transpose();
		fu[i] = fu_arg[i];
		Fu_bar[i] = Fu[i]*K;		// the first one is NOT used
		Fu_bar_transp[i] = Fu_bar[i].transpose();
	}
	
	F_xTheta = F_xTheta_arg;
	F_xTheta_transp = F_xTheta_arg.transpose();
	F_theta = F_theta_arg;
	F_theta_transp = F_theta_arg.transpose();
	f_xTheta = f_xTheta_arg;
	num_constr = _N*(_nInp + _nSt) + _nF_xTheta;
	
	s = s_arg;
	
	// do some preliminary Matrix calculations	
	S = K.transpose() * R;
	S_transp = S.transpose();
	Q_bar = S * K;
	A_bar = A + B*K;
	A_bar_transp = A_bar.transpose();	
	
	// vvv build the matrices Fx_hat(_transp), Fu_hat(_transp), F_xTheta_hat(_transp), F_theta_hat(_transp), Fu_bar(_transp) for PhaseI
	Matrix<Type, _nInp, 1> ones_nInp;
	ones_nInp.setConstant(-1);
	Matrix<Type, _nSt, 1> ones_nSt;
	ones_nSt.setConstant(-1);
	
	Fu_hat[0] << Fu[0], ones_nInp;
	Fu_hat_transp[0] = Fu_hat[0].transpose();
	for (int i=1; i<=_N-1; i++)
	{
		Fx_hat[i-1] << Fx[i-1], ones_nSt;
		Fx_hat_transp[i-1] = Fx_hat[i-1].transpose();
		
		Fu_bar_hat[i-1].setZero();
		Fu_bar_hat[i-1].block(0,0,_nInp,_n) = Fu_bar[i];
		Fu_bar_hat_transp[i-1] = Fu_bar_hat[i-1].transpose();
		
		Fu_hat[i] << Fu[i], ones_nInp;
		Fu_hat_transp[i] = Fu_hat[i].transpose();
	}
	Fx_hat[_N-1] << Fx[_N-1], ones_nSt;
	Fx_hat_transp[_N-1] = Fx_hat[_N-1].transpose();
	
	F_xTheta_hat.setZero();
	F_xTheta_hat.block(0,0,_nF_xTheta,_n) = F_xTheta;
	F_xTheta_hat_transp = F_xTheta_hat.transpose();
	
	F_theta_hat << F_theta, Matrix<Type, _nF_xTheta, 1>::Constant(_nF_xTheta, 1, -1);
	F_theta_hat_transp = F_theta_hat.transpose();
	// ^^^ built the matrices
	
	z_warm.setConstant(100);
	
}


// step function returns optimal input
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
Matrix<Type, _m, 1> LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega>::step(const Matrix<Type, _n, _n> &Lm_arg, const Matrix<Type, _n, _m> &Mm_arg, const Matrix<Type, _n, 1> &tm_arg,
																									   const Matrix<Type, _n, 1> &x_hat_arg,
																									  const Matrix<Type, _n, 1> x_star_arg[]
												   													  )
{
	// initialization
	Lm = Lm_arg;
	Mm = Mm_arg;
	tm = tm_arg;
	x_hat = x_hat_arg;
	kappa = kappa_start;
	x_star = x_star_arg;
	tm_tilde = tm_arg + s;
	Am_tilde = A + Lm;
	Am_tilde_transp = Am_tilde.transpose();
	Bm = B + Mm;
	Bm_transp = Bm.transpose();
	Bm_bar = Bm * K;
	Bm_bar_transp = Bm_bar.transpose();
	
	//update z_warm from previous z_warm
	z_warm.template segment<(_N-1)*(_m+_n+_n)>(0) = z_warm.template segment<(_N-1)*(_m+_n+_n)>(offset);
	
	// ONLY TEMPORARILY
	// z_warm << -2 ,    3  ,   3  ,   1 ,   -2 ,    2 ,   -3 ,    0  ,  -1  ,   0  ,   5 ,    3  ,   5  ,  -3  ,   0  ,  -4  ,   3  ,   1  ,   4  ,   5  ,   4 ,   -1,
		// -5  ,   0  ,  -3 ,   -3  ,  -2  ,  -4  ,   2  ,   2   ,  0  ,  -2  ,   3  ,   1 ,    5  ,   4  ,  -1  ,   0   , -2   ,  1  ,   3  ,   2 ,   -4   ,  3,
		// -5   , -1   ,  2   ,  3  ,  -1   ,  2;
	z_warm.setZero();
	
	compRQ();	// compute u_star, x_star -> cost matrices 
	z = z_warm;
	
	nu.setConstant(1);	// any nu good to initialize
	
	// tests, if given w_warm satisfies (P*z_warm<h), otherwise a suitable z_warm is computed using phaseI()
	// testPos(): returns 1 if inequality is satisfied
	//					  0 if inequality not satisfied
	if( !testPos() )	// testPos() operates on z
	{
			if (!phaseI())	// if phaseI() == 0 => Phase I was not successful
			{
				cerr << "Phase I was unsuccessful" << endl;
				cerr << "process aborted." << endl;
				return K*x_hat + z.template segment<_m>(0);
			}		
	}
	else	// sets new z_warm and an appropriate nu
	{
		cout << "chosen z_warm is in domain." << endl;
	}
	// return K*x_hat + z.template segment<_m>(0);

	kappa = kappa_start;
	z = z_warm;		// the (new) z_warm satisfies inequality
	compD();	// compute d-vector
	
	Matrix<Type, _N*(_m + _n + _n) + _m, 1> z_orig;	// stores the "old z", z now stores z = z_orig + t*dz
	Matrix<Type, 2*_N*_n, 1> nu_orig;	// ditto
	double resNormSq_orig;	// squared norm of original resudia
	double resNormSq;	// squared norm of residua
	double resNormSqRp;	// squared norm of primal residuum
	bool cond;
	int i;		// counter for backtracking line search
	int j;		// number of newton iterations for a fixed kappa
	int loopCounter = 0;
	int i_kappa = 0;	// delete it later
	do 			// kappa changes
	{
		i_kappa++;	
		// cout << "============= kappa: " << kappa << ", i_kappa: " << i_kappa << " ================" << endl;
		compRdRp();		// z is already in domain, using new kappa to compute primal and dual residua
		j = 0;
		do 		// kappa fixed
		{
			j++;	// number of Newton steps for fixed kappa
			compDzDnu();
			
			// --------- line search part --------------
			// get smaller and smaller t until following two conditions are satisfied
			// 1) ||r_d(z+t*dz,nu+t*dnu)||_2 + ||r_p(z+t*dz,nu+t*dnu)||_2 <= (1-alpha*t)*(||r_d||_2 + ||r_p||_2)
			// 2) P*(z+t*dz) < h
			t = compT();	// retrieve largest t s.t. P*(z+t*dz) < h			
			z_orig = z;	// stores the "old z", z stores z = z_tmp + t*dz
			nu_orig = nu;
			
			resNormSq_orig = r_d.squaredNorm() + r_p.squaredNorm();
			cond = 1;
			i = 0;
			loopCounter++;
			
			while (cond)	// backtracking line search
			{
				i++;
				z = z_orig + t*dz;
				nu = nu_orig + t*dnu;
				compD();
				compRdRp();		// class variables r_d and r_p are updated
				resNormSqRp = r_p.squaredNorm();
				resNormSq = r_d.squaredNorm() + resNormSqRp;
				cond = ( (resNormSq > (1-alpha_ls*t)*(1-alpha_ls*t)*resNormSq_orig) );
					
				t = beta_ls*t;
				if (t <= eps_ls)
					break;	// stop, no improvement assumed
			}	// z and nu have automatically been updated
			// cout << "line search needed " << i << " steps with step size t: " << t/beta_ls << endl << endl;
			
			if (j >= n_iter && i_kappa > 0)		// >0 includes first step i_kappa==1, even though it might need more steps
			{
				// cout << "more than " << n_iter <<" steps required" << endl;
				// cout << "loopCounter: " << loopCounter << endl;
				break;
			}
		} while( (resNormSq > eps_ntSq) || (resNormSqRp > eps_normRpSq) );
		kappa = kappa*mu;
	} while( kappa/mu*num_constr > eps_barrier );
	
	cout << " =====> computed optimal z_vector:" << endl << setprecision (15) << z << endl << endl;
	cout << setprecision (15) << "kappa: " << kappa << " after i_kappa: " << i_kappa << endl;
	cout << setprecision (15) << "eps_barrier: " << eps_barrier << endl;
	cout << setprecision (15) << "kappa/mu*num_constr: " << kappa/mu*num_constr << endl << endl;
	cout << setprecision (15) << "loopCounter: " << loopCounter << endl << endl;
	
	z_warm = z;
	return K*x_hat + z.template segment<_m>(0);
}


// ------------ function tests if z_warm satisfies P*z < h --------
// ------- might consider computing entire matrix check and then iterate for neg. elements ----
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: testPos()
{	
	// special treatment at beginning
	c_tmp = z.template segment<_m>(0);
	checkU = (fu[0] - Fu_bar[0]*x_hat) - (Fu[0]*c_tmp);	// should be >0
	
	for (int j=1; j <= _nInp; j++)
	{
		if (checkU[j-1] <= 0)
			return 0;
	}
	
	// general treatment in the middle, class variable: offset = _n + _n + _m;
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
				return 0;
		}
		
		checkU = fu[i] - (Fu_bar[i]*x_bar_tmp + Fu[i]*c_tmp);
		for (int j=1; j<=_nInp; j++)
		{
			if (checkU[j-1] <= 0)
				return 0;
		}
	}
	
	// special case for last blocks
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);	// depends on where it is
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	checkX = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	for (int j=1; j<= _nSt; j++)
	{
		if (checkX[j-1] <= 0)
			return 0;
	}
	
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);	// depends on position of invariant set
	checkTheta = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if (checkTheta[j-1] <= 0)
			return 0;
	}
	
	return 1;	// no problems detected
}


// ------------ function finds z_warm such that P*z_warm < h -----------
// ------------ similar to PhaseII, heuristic method ------------
// solve the problem: min{gamma_0 + gamma_2N}, s.t. {P*z-h < gamma, C*z = b, gamma_0=gamma_1=...=gamma_2N-1}
// reformulate as: min{z_hat*H_hat*z + g_hat'*z_hat}, s.t. {P_hat*z_hat < h, C_hat * z_hat = b_hat}
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: phaseI()
{
	// cout << "---------- Pz_warm < h not satisfied, starting phaseI() -----------" << endl;
	nu_hat.setConstant(1);	// any nu good to initialize
	
	// for test purposes, n=5, m=2, N = 10
	 nu_hat << 7, 5, 7, 10, 5, 3, 8, 2, 2, 5, 6, 6, 1, 1, 3, 7, 4, 8, 9, 4, 1, 6, 7, 3, 5, 2, 9, 5, 5, 2, 1, 1, 8, 8, 2, 9, 7, 7, 10, 0, 4, 5, 6, 1,
	 7, 1, 6, 7, 5, 6, 2, 6, 4, 1, 4, 3, 2, 6, 5, 4, 2, 9, 9, 3, 6, 8, 10, 9, 4, 2, 8, 1, 2, 4, 8, 6, 6, 9, 6, 0, 8, 6, 5, 3, 3, 5, 2, 0,
	 2, 3, 2, 6, 8, 8, 8, 7, 9, 8, 2, 10, 1, 7, 2, 2, 3, 2, 5, 10, 4, 6, 6, 7, 4, 9, 9, 8, 3, 8, 1;
	
	
	compMatrices_PhaseI();	// computes Fx_hat, Fu_hat, Fu_bar_hat, g_hat_c, g_hat_x, Bm_hat, Identity_hat...
	difference = 10;	// h-P*z_warm = difference
	compGamma_PhaseI();	// computes gamma
	
	// vvv build z_hat, offset2 = _n+_n+_m+2;
	for (int i=1; i<= _N; i++)
	{
		z_hat.template segment<_m>( (i-1)*offset2 ) = z.template segment<_m>( (i-1)*offset );
		z_hat.template segment<1>((i-1)*offset2+_m) = gamma.template segment<1>(0);	// same
		z_hat.template segment<_n+_n>( (i-1)*offset2+_m+1) = z.template segment<_n+_n>( (i-1)*offset+_m );
		z_hat.template segment<1>((i-1)*offset2+_m+1+_n+_n) = gamma.template segment<1>(0);
	}	
	z_hat.template segment<_m>(_N*offset2) = z.template segment<_m>(_N*offset);
	z_hat.template segment<1>(_N*offset2+_m) = gamma.template segment<1>(2*_N);
	// cout << setprecision(30) << "z_hat:" << endl << z_hat << endl << endl;
	// ^^^ computed z_hat
		
	compD_hat_PhaseI();	// stored in d 
	double resNormSq_orig;
	double resNormSq;
	double resNormSqRp;
	int j;
	int loopCounter = 0;
	int i_kappa = 0;
	kappa = kappa_start_PhaseI;
	bool cond;
	int i;

	do
	{
		i_kappa++;
		// cout << "============= kappa PhaseI: " << kappa << ", i_kappa: " << i_kappa << " ================" << endl;
		compRdRp_PhaseI();		// computes r_p and r_d_hat
		// cout << setprecision(30) << "r_p_hat" << endl << r_p_hat << endl << endl;
		// cout << setprecision(30) << "r_d_hat" << endl << r_d_hat << endl << endl;
		// return 0;
		j = 0;
		do  	// kappa fixed 
		{
			j++;
			// cout << "----------- kappa fixed PhaseI(), round: " << j << "--------------" << endl << endl;			
			compDz_hatDnu_hat_PhaseI();	// compute dz_hat and dnu
			// cout << setprecision(30) << "dz_hat" << endl << dz_hat << endl << endl;
			// cout << "dnu_hat" << endl << dnu_hat << endl << endl;
			// return 0;
			// ---- line search part --------
			// decrease t until both conditions are satisfied
			// 1) ||res(z+t*dz,nu+t*dnu)||_2 <= (1-alpha*t)*(||r_d||_2 + ||r_p||_2)
			// 2) P*(z+t*dz) < h
			t = compT_PhaseI();		// computes largest t: P_hat*(z_hat+t*dz_hat) < h
			// cout << "t: " << endl << t << endl << endl;
			z_hat_orig = z_hat;
			nu_hat_orig = nu_hat;
			resNormSq_orig = r_d_hat.squaredNorm() + r_p_hat.squaredNorm();

			cond = 1;
			i = 0;
			loopCounter++;
			while (cond)
			{
				i++;
				z_hat = z_hat_orig + t*dz_hat;	// guaranteed to be feasible, b/c t is small enough
				// stop if gamma < 0, z_warm is computed automatically
				if( testNeg_PhaseI() )
				{
					// cout << setprecision(30) << "new z_warm is: " << endl << z_warm << endl << endl;
					cout << "PhaseI needed " << loopCounter << " Newton steps to find z_warm." << endl << endl;
					return 1;
				}
				nu_hat = nu_hat_orig + t*dnu_hat;
				compD_hat_PhaseI();	// stored in d
				compRdRp_PhaseI();
				resNormSqRp = r_p_hat.squaredNorm();
				resNormSq = r_d_hat.squaredNorm() + resNormSqRp;
				cond = ( (resNormSq > (1-alpha_ls*t)*(1-alpha_ls*t)*resNormSq_orig) );
				t = beta_ls*t;
			}	
			if (j >= n_iter_PhaseI)
			{
				// cout << " PhaseI more than " << n_iter_PhaseI << " steps for kappa = " << kappa << endl;
				// cout << "loopCounter: " << loopCounter << endl;
				break;
			}
		} while( (resNormSq > eps_ntSq) || (resNormSqRp > eps_normRpSq) );
		kappa = kappa*mu;
	} while( kappa/mu*num_constr > eps_barrier );
	cerr << "***** phaseI() did NOT converge." << endl;
	return 0;
}


// ------------ function computes primal and dual residua  -------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compRdRp()
{	// ---------- compute primal dual r_p = C*z - b, done two blocks
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
	
	// ------------- compute dual r_d = 2*H*z + g + kappa*P'*d + C'*nu;
	// handle first case separately
	r_d.template segment<_m>(0) = 2*R*z.template segment<_m>(0) + (r_vec[0]+2*S_transp*x_hat) + 
				kappa*Fu_transp[0]*d.template segment<_nInp>(0) - Bm_transp*nu.template segment<_n>(0) - B_transp*nu.template segment<_n>(_n);
						
	// handle the cases in the middle, without the end, three block to deal with in each round
	int offset1 = 2*_n;	// offset required in nu for C'*nu
	for (int i=1; i<= _N-1; i++)
	{
		r_d.template segment<_n>(_m+(i-1)*offset) = 2*Q_tilde*z.template segment<_n>(_m+(i-1)*offset) + q_tilde_vec[i-1] +
				nu.template segment<_n>((i-1)*offset1) - Am_tilde_transp*nu.template segment<_n>((i-1)*offset1+_n+_n);
		
		if (i != _pos_omega)	
		{	
			r_d.template segment<_n>(_m+(i-1)*offset+_n) = 2*Q_bar*z.template segment<_n>(_m+(i-1)*offset+_n) + 2*S*z.template segment<_m>(_m+(i-1)*offset+_n+_n) + q_bar_vec[i] +
					kappa*Fx_transp[i-1]*d.template segment<_nSt>(_nInp+(i-1)*(_nInp+_nSt)) + kappa*Fu_bar_transp[i]*d.template segment<_nInp>(i*(_nInp+_nSt)) + 
					nu.template segment<_n>((i-1)*offset1+_n) - Bm_bar_transp*nu.template segment<_n>((i-1)*offset1+_n+_n) - A_bar_transp*nu.template segment<_n>((i-1)*offset1+_n+_n+_n);
		}
		else	// must add the additional term: F_xTheta'*d(.)
		{
			r_d.template segment<_n>(_m+(_pos_omega-1)*offset+_n) = 2*Q_bar*z.template segment<_n>(_m+(_pos_omega-1)*offset+_n) + 2*S*z.template segment<_m>(_m+(_pos_omega-1)*offset+_n+_n) + q_bar_vec[i] +
					kappa*Fx_transp[_pos_omega-1]*d.template segment<_nSt>(_nInp+(i-1)*(_nInp+_nSt)) + kappa*Fu_bar_transp[_pos_omega]*d.template segment<_nInp>(i*(_nInp+_nSt)) + 
					nu.template segment<_n>((_pos_omega-1)*offset1+_n) - Bm_bar_transp*nu.template segment<_n>((_pos_omega-1)*offset1+_n+_n) - A_bar_transp*nu.template segment<_n>((_pos_omega-1)*offset1+_n+_n+_n)
					+ kappa*F_xTheta_transp * d.template segment<_nF_xTheta>(num_constr-_nF_xTheta);
		}
		
		r_d.template segment<_m>(_m+(i-1)*offset+_n+_n) = 2*S_transp*z.template segment<_n>(_m+(i-1)*offset+_n) + 2*R*z.template segment<_m>(_m+(i-1)*offset+_n+_n) + r_vec[i] +
				// kappa*Fu_transp[i]*d.segment(tmp,Fu_rows[i]) - 
				kappa*Fu_transp[i]*d.template segment<_nInp>(i*(_nInp+_nSt)) - 
				Bm_transp*nu.template segment<_n>((i-1)*offset1+_n+_n) - B_transp*nu.template segment<_n>((i-1)*offset1+_n+_n+_n);		
	}
	r_d.template segment<_n>(_m+(_N-1)*offset) = 2*Q_tilde_f*z.template segment<_n>(_m+(_N-1)*offset) + q_tilde_vec[_N-1] + nu.template segment<_n>((_N-1)*offset1);
	
	if(_pos_omega == _N)
	{
		r_d.template segment<_n>(_m+(_N-1)*offset+_n) = kappa*Fx_transp[_N-1]*d.template segment<_nSt>(_nInp+(_N-1)*(_nInp+_nSt)) + kappa*F_xTheta_transp*d.template segment<_nF_xTheta>(_N*(_nInp+_nSt)) + 
				nu.template segment<_n>((_N-1)*offset1+_n);
	}
	else	//standard
	{
		r_d.template segment<_n>(_m+(_N-1)*offset+_n) = kappa*Fx_transp[_N-1]*d.template segment<_nSt>(_nInp+(_N-1)*(_nInp+_nSt)) + nu.template segment<_n>((_N-1)*offset1+_n);
	}
	r_d.template segment<_m>(_m+(_N-1)*offset+_n+_n) = kappa * F_theta_transp*d.template segment<_nF_xTheta>(_N*(_nSt+_nInp));
}


// ------------ function computes vector (d)_i=(h-P*z)_i  --------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compD()
{
	// special treatment at beginning
	c_tmp = z.template segment<_m>(0);
	checkU = (fu[0] - Fu[0]*K*x_hat) - (Fu[0]*c_tmp);	// should be >0
	for (int j=1; j <= _nInp; j++)
	{
		d[j-1] = 1/checkU[j-1];
	}

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
		
		checkU = fu[i] - (Fu[i]*K*x_bar_tmp + Fu[i]*c_tmp);
		for (int j=1; j<=_nInp; j++)
		{
			d[_nInp+(i-1)*(_nInp+_nSt)+_nSt+j-1] = 1/checkU[j-1];
		}
	}
	
	// special case for last blocks
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	
	checkX = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	for (int j=1; j<= _nSt; j++)
	{
		d[_nInp+(_N-1)*(_nSt+_nInp)+j-1] = 1/checkX[j-1];
	}
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);
	checkTheta = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		d[_N*(_nSt+_nInp)+j-1] = 1/checkTheta[j-1];
	}
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
	/*
	DiagonalMatrix<Type, Dynamic> d_diag[2*_N+1];
	for (int i=0 ; i <= _N-1 ; i++)
	{
		d_diag[2*i].diagonal() = d.template segment<_nInp>(i*(_nInp + _nSt));
		d_diag[2*i+1].diagonal() = d.template segment<_nSt>(_nInp+i*(_nInp + _nSt));
	}
	d_diag[2*_N].diagonal() = d.template segment<_nF_xTheta>(_N*(_nInp + _nSt));
	*/
	DiagonalMatrix<Type, _nInp> d_diagU[_N];
	DiagonalMatrix<Type, _nSt> d_diagX[_N];
	DiagonalMatrix<Type, _nF_xTheta> d_diagTheta;
	for (int i=0 ; i <= _N-1 ; i++)
	{
		d_diagU[i].diagonal() = d.template segment<_nInp>(i*(_nInp + _nSt));
		d_diagX[i].diagonal() = d.template segment<_nSt>(_nInp+i*(_nInp + _nSt));
	}
	d_diagTheta.diagonal() = d.template segment<_nF_xTheta>(_N*(_nInp + _nSt));
	
	// --------------- compute elements of Phi, i.e. elements of Omicron, Pi = 2*Q_tilde, Rho, Sigma		
	Matrix<Type, _m, _m> eyeM;
	eyeM.setIdentity();
	Matrix<Type, _n, _n> eyeN;
	eyeN.setIdentity();
	
	// special treatment at the beginning
	// Omicron[0] = 2*R + kappa*Fu_transp[0]*d_diag[0]*d_diag[0]*Fu[0]	+ reg * eyeM;
	Omicron[0] = 2*R + kappa*Fu_transp[0]*d_diagU[0]*d_diagU[0]*Fu[0]	+ reg * eyeM;
		
	// do the rest by computing three block and three block
	for (int i=1; i <= _N-1; i++)
	{
		if (i != _pos_omega)
		{
			// Rho[i-1] = 2*Q_bar + kappa * ( Fx_transp[i-1]*d_diag[2*(i-1)+1]*d_diag[2*(i-1)+1]*Fx[i-1]  + Fu_bar_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu_bar[i] ) + reg * eyeN;
			Rho[i-1] = 2*Q_bar + kappa * ( Fx_transp[i-1]*d_diagX[i-1]*d_diagX[i-1]*Fx[i-1]  + Fu_bar_transp[i]*d_diagU[i]*d_diagU[i]*Fu_bar[i] ) + reg * eyeN;
		}
		else	// i == _pos_omega
		{
			// Rho[_pos_omega-1] = 2*Q_bar + kappa * ( Fx_transp[_pos_omega-1]*d_diag[2*(_pos_omega-1)+1]*d_diag[2*(_pos_omega-1)+1]*Fx[i-1]  + Fu_bar_transp[_pos_omega]*d_diag[2*_pos_omega]*d_diag[2*_pos_omega]*Fu_bar[_pos_omega] 
				// + F_xTheta_transp * d_diag[2*_N]*d_diag[2*_N]*F_xTheta)	+ reg * eyeN;
			Rho[_pos_omega-1] = 2*Q_bar + kappa * ( Fx_transp[_pos_omega-1]*d_diagX[_pos_omega-1]*d_diagX[_pos_omega-1]*Fx[i-1]  + Fu_bar_transp[_pos_omega]*d_diagU[_pos_omega]*d_diagU[_pos_omega]*Fu_bar[_pos_omega] 
				+ F_xTheta_transp * d_diagTheta*d_diagTheta*F_xTheta)	+ reg * eyeN;
		}
		// Sigma[i-1] = 2*S + kappa * ( Fu_bar_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu[i] );
		Sigma[i-1] = 2*S + kappa * ( Fu_bar_transp[i]*d_diagU[i]*d_diagU[i]*Fu[i] );
		// Omicron[i] = 2*R + kappa * ( Fu_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu[i] ) + reg * eyeM;
		Omicron[i] = 2*R + kappa * ( Fu_transp[i]*d_diagU[i]*d_diagU[i]*Fu[i] ) + reg * eyeM;
	}
	
	// special treatment for last block
	if (_pos_omega == _N)
	{
		// Rho[_N-1] = kappa * ( Fx_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx[_N-1] + F_xTheta_transp*d_diag[2*_N]*d_diag[2*_N]*F_xTheta ) + reg * eyeN;
		Rho[_N-1] = kappa * ( Fx_transp[_N-1]*d_diagX[_N-1]*d_diagX[_N-1]*Fx[_N-1] + F_xTheta_transp*d_diagTheta*d_diagTheta*F_xTheta ) + reg * eyeN;
	}
	else	// considered in loop above
	{
		// Rho[_N-1] = kappa * ( Fx_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx[_N-1]) + reg * eyeN;
		Rho[_N-1] = kappa * ( Fx_transp[_N-1]*d_diagX[_N-1]*d_diagX[_N-1]*Fx[_N-1]) + reg * eyeN;
	}
	// Sigma[_N-1] = kappa * ( F_xTheta_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta );	// independent of _pos_omega, represents the off-diag matrix
	Sigma[_N-1] = kappa * ( F_xTheta_transp*d_diagTheta*d_diagTheta*F_theta );	// independent of _pos_omega, represents the off-diag matrix
	// Omicron[_N] = kappa * ( F_theta_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta ) + reg * eyeM;
	Omicron[_N] = kappa * ( F_theta_transp*d_diagTheta*d_diagTheta*F_theta ) + reg * eyeM;
	
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
	
	// decompose Phi = L_Phi*L_Phi'
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
		U[5+(_N-2)*3] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( /*zero*/ - LSigma_offDiag[_N-1]*U_bar[6+(_N-2)*5] );
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
	
	// special treatment for L22, L32, L42
	L_diag[1].compute( Y[1][1]-L_offDiag[0][0]*L_offDiag_transp[0][0] );
	L_diag_transp[1] = L_diag[1].matrixLLT().transpose();
	
	L_offDiag_transp[0][1] = L_diag[1].matrixLLT().triangularView<Lower>().solve(Y[1][2]-L_offDiag[0][0]*L_offDiag_transp[1][0]);
	L_offDiag[0][1] = L_offDiag_transp[0][1].transpose();
	
	L_offDiag_transp[1][1] = L_diag[1].matrixLLT().triangularView<Lower>().solve(Y[0][3]);
	L_offDiag[1][1] = L_offDiag_transp[1][1].transpose();
	
	// cases in the middle
	for (int i = 1; i <= 2*_N-4; i++)
	{
		L_diag[i+1].compute( Y[2][i+1] - L_offDiag[1][i-1]*L_offDiag_transp[1][i-1] - L_offDiag[0][i]*L_offDiag_transp[0][i] );
		L_diag_transp[i+1] = L_diag[i+1].matrixLLT().transpose();
		L_offDiag_transp[0][i+1] = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( Y[1][i+2] - L_offDiag[0][i]*L_offDiag_transp[1][i] );
		L_offDiag[0][i+1] = L_offDiag_transp[0][i+1].transpose();
		
		L_offDiag_transp[1][i+1] = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( Y[0][i+3] );
		L_offDiag[1][i+1] = L_offDiag_transp[1][i+1].transpose();
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
		delta.template segment<_n>(_n+i*_n) = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( -beta.template segment<_n>(_n+i*_n) - L_offDiag[1][i-1]*delta.template segment<_n>((i-1)*_n) - L_offDiag[0][i]*delta.template segment<_n>(i*_n) );
	
	// 2) now, solve for L'*Dnu = delta
	dnu.template segment<_n>(2*_n*_N - _n) = L_diag_transp[2*_N-1].template triangularView<Upper>().solve(delta.template segment<_n>(2*_n*_N - _n) );
	dnu.template segment<_n>(2*_n*_N - _n - _n) = L_diag_transp[2*_N-2].template triangularView<Upper>().solve( delta.template segment<_n>(2*_n*_N - _n - _n) - L_offDiag_transp[0][2*_N-2]*dnu.template segment<_n>(2*_n*_N - _n) );
	
	//remaining cases are regular
	for (int i=1; i<=2*_N-2; i++)
		dnu.template segment<_n>(2*_n*_N-(i+2)*_n) = L_diag_transp[2*_N-(i+2)].template triangularView<Upper>().solve( delta.template segment<_n>(2*_n*_N-(i+2)*_n) - L_offDiag_transp[0][2*_N-(i+2)]*dnu.template segment<_n>(2*_n*_N-(i+1)*_n) - L_offDiag_transp[1][2*_N-(i+2)]*dnu.template segment<_n>(2*_n*_N-i*_n)  );
	//cout << "dnu" << endl << dnu << endl << endl;
}


// ------------ function computes Phi * dz = -r_d - C' * dnu ---------------------
// ------------ implementing --------------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDz()
{	
	// computed in two parts
	// 1. tmp = -r_d - C' * dnu
	// 2. L_Phi*L_Phi'*dz = tmp
	// 3. L_Phi*tmp1 = tmp
	// 4. L_Phi'*dz = tmp1
	
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
	
	
	// 3. L_Phi*tmp1 = tmp
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
	
	
	// 4. L_Phi'*dz = tmp1
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
	
	double t = 1;	// stores smallest variable
	double t_tmp;	// t_tmp = min(t_vec)
	
	c_tmp = z.template segment<_m>(0);
	dc_tmp = dz.template segment<_m>(0);
	checkU = (fu[0] - Fu_bar[0]*x_hat) - (Fu[0]*c_tmp);	// should be >0
	dcheckU = Fu[0] * dc_tmp;
	t_vecU.setConstant(1);
	for (int j=1; j <= _nInp; j++)
	{
		if (dcheckU[j-1] > 0)	// neg. cases not interesting
			t_vecU[j-1] = checkU[j-1]/dcheckU[j-1];
	}
	t_tmp = t_vecU.minCoeff();
	if (t_tmp < t)
		t = t_tmp;
	
	
	// class variable: offset = _n + _n + _m;
	// general treatment in the middle 
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	// compute (h - P*z)_i
		x_bar_tmp = z.template segment<_n>((i-1)*offset+_m+_n);
		dx_bar_tmp = dz.template segment<_n>((i-1)*offset+_m+_n);
		c_tmp = z.template segment<_m>((i-1)*offset+_m+_n+_n);
		dc_tmp = dz.template segment<_m>((i-1)*offset+_m+_n+_n);
		checkX = fx[i-1] - Fx[i-1]*x_bar_tmp;
		dcheckX = Fx[i-1]*dx_bar_tmp;
		t_vecX.setConstant(1);
		for (int j=1; j<= _nSt; j++)
		{
			if(dcheckX[j-1]>0)
				t_vecX[j-1] = checkX[j-1]/dcheckX[j-1];
		}
		t_tmp = t_vecX.minCoeff();
		if (t_tmp < t)
			t = t_tmp;
		
		checkU = fu[i] - (Fu_bar[i]*x_bar_tmp + Fu[i]*c_tmp);
		dcheckU = Fu_bar[i]*dx_bar_tmp + Fu[i]*dc_tmp;
		t_vecU.setConstant(1);
		for (int j=1; j<=_nInp; j++)
		{
			if(dcheckU[j-1]>0)
				t_vecU[j-1] = checkU[j-1]/dcheckU[j-1];
		}
		t_tmp = t_vecU.minCoeff();
		if (t_tmp < t)
			t = t_tmp;
	}
	
	// special case for last blocks
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);
	dx_bar_tmp = dz.template segment<_n>((_N-1)*offset+_m+_n);
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	dc_tmp = dz.template segment<_m>((_N-1)*offset+_m+_n+_n);
	checkX = fx[_N-1] - (Fx[_N-1]*x_bar_tmp);
	dcheckX = Fx[_N-1]*dx_bar_tmp;
	t_vecX.setConstant(1);
	for (int j=1; j<= _nSt; j++)
	{
		if (dcheckX[j-1]>0)
			t_vecX[j-1] = checkX[j-1]/dcheckX[j-1];
	}
	t_tmp = t_vecX.minCoeff();
	if (t_tmp < t)
		t = t_tmp;
	
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);
	dx_bar_tmp = dz.template segment<_n>((_pos_omega-1)*offset+_m+_n);
	checkTheta = f_xTheta - (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	dcheckTheta = F_xTheta*dx_bar_tmp + F_theta*dc_tmp;
	t_vecTheta.setConstant(1);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if(dcheckTheta[j-1]>0)
			t_vecTheta[j-1] = checkTheta[j-1]/dcheckTheta[j-1];
	}
	t_tmp = t_vecTheta.minCoeff();
	if (t_tmp < t)
		t = t_tmp;
	
	// return the result
	if (t == 1)
		return 1;
	else
		return 0.99*t;	// guarantees strict feasibility
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
		u_star[i] = svd.solve(x_star[i] - Am_tilde*x_star[i-1] - tm_tilde);
	
	/*
	for (int i = 0; i <= _N-1; i++)
		cout << "u_star[" << i << "]" << endl << u_star[i] << endl << endl;
	*/
	
	// compute the vectors q_bar_vec[]; q_tilde_vec[]; r_vec[]
	q_tilde_vec[0] = -2*Q_tilde*x_star[0];
	r_vec[0] = -2*R*u_star[0];
	// q_bar_vec[0] is never used
	for (int i = 1; i <= _N-2; i++)		// be careful how to use it
	{
		q_tilde_vec[i] = -2*Q_tilde*x_star[i];
		r_vec[i] = -2*R*u_star[i];
		q_bar_vec[i] = K_transp*r_vec[i];	
	}
	q_tilde_vec[_N-1] = -2*Q_tilde_f*x_star[_N-1];
	r_vec[_N-1] = -2*R*u_star[_N-1];
	q_bar_vec[_N-1] = K_transp*r_vec[_N-1];
	
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
	
	// vvv definitions of Bm_hat, B_hat, B_bar_hat, A_bar_hat, Identity_hat
	Bm_hat.setZero();
	Bm_hat.template block<_n,_m>(0,0) = Bm;
	Bm_hat.bottomRightCorner(1,1).setConstant(1);
	Bm_hat_transp = Bm_hat.transpose();
	
	Bm_hat0.setZero();
	Bm_hat0.template block<_n,_m>(0,0) = Bm;
	Bm_hat0_transp = Bm_hat0.transpose();
	
	B_hat.setZero();
	B_hat.template block<_n,_m>(0,0) = B;
	B_hat.bottomRightCorner(1,1).setConstant(-1);
	B_hat_transp = B_hat.transpose();

	Identity1_hat.setZero();
	Identity1_hat.topLeftCorner(_n,_n).setIdentity();	// or: Identity_hat[i-1].topLeftCorner(_n,_n).setIdentity();
	Identity1_hat_transp = Identity1_hat.transpose();
	
	Identity2_hat.setZero();
	Identity2_hat.topLeftCorner(_n,_n).setIdentity();
	Identity2_hat.bottomRightCorner(1,1).setConstant(-1);
	Identity2_hat_transp = Identity2_hat.transpose();
	
	Bm_bar_hat.setZero();
	Bm_bar_hat.template block<_n,_n>(0,0) = Bm_bar;
	Bm_bar_hat.bottomRightCorner(1,1).setConstant(-1);
	Bm_bar_hat_transp = Bm_bar_hat.transpose();
	
	A_bar_hat.setZero();
	A_bar_hat.template block<_n,_n>(0,0) = A_bar;
	A_bar_hat_transp = A_bar_hat.transpose();
	
	Am_tilde_hat.setZero();
	Am_tilde_hat.template block<_n,_n>(0,0) = Am_tilde;
	Am_tilde_hat_transp = Am_tilde_hat.transpose();
	
	
	/*
	cout << "Bm_hat" << endl << Bm_hat << endl << endl;
	cout << "Bm_hat_transp" << endl << Bm_hat_transp << endl << endl;
	cout << "B_hat" << endl << B_hat << endl << endl;
	cout << "B_hat_transp" << endl << B_hat_transp << endl << endl;
	cout << "Identity1_hat" << endl << Identity1_hat << endl << endl;
	cout << "Identity1_hat_transp" << endl << Identity1_hat_transp << endl << endl;
	cout << "Identity2_hat" << endl << Identity2_hat << endl << endl;
	cout << "Identity2_hat" << endl << Identity2_hat << endl << endl;
	cout << "Bm_bar_hat" << endl << Bm_bar_hat << endl << endl;
	cout << "Bm_bar_hat_transp" << endl << Bm_bar_hat_transp << endl << endl;
	cout << "A_bar_hat" << endl << A_bar_hat << endl << endl;
	cout << "A_bar_hat_transp" << endl << A_bar_hat_transp << endl << endl;
	cout << "Am_tilde_hat" << endl << Am_tilde_hat << endl << endl;
	cout << "Am_tilde_hat_transp" << endl << Am_tilde_hat_transp << endl << endl;
	*/
}


// ------------ function computes a feasibe gamma needed by phaseI ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compGamma_PhaseI()
{
	double gamma0;	// stores the largest gamma necessary for the state and input constraints, i.e. s_0, ... , s_2N-1
	double gamma_tmp;		// variable used to comporarily compute variable gamma0
	
	// special treatment at beginning	// c_tmp = z.segment(0,_m);
	c_tmp = z.template segment<_m>(0);
	checkU = -(fu[0] - Fu[0]*K*x_hat) + (Fu[0]*c_tmp);	// should be >0
	gamma0 = checkU.maxCoeff();
	
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		x_bar_tmp = z.template segment<_n>((i-1)*offset+_m+_n);
		c_tmp = z.template segment<_m>((i-1)*offset+_m+_n+_n);		
		checkX = -fx[i-1] + Fx[i-1]*x_bar_tmp;
		gamma_tmp = checkX.maxCoeff();
		if (gamma0 < gamma_tmp)
			gamma0 = gamma_tmp;
			
		checkU = -fu[i] + (Fu[i]*K*x_bar_tmp + Fu[i]*c_tmp);
		gamma_tmp = checkU.maxCoeff();
		if (gamma0 < gamma_tmp)
			gamma0 = gamma_tmp;
	}
	x_bar_tmp = z.template segment<_n>((_N-1)*offset+_m+_n);
	c_tmp = z.template segment<_m>((_N-1)*offset+_m+_n+_n);
	checkX = -fx[_N-1] + (Fx[_N-1]*x_bar_tmp);
	gamma_tmp = checkX.maxCoeff();
	if (gamma0 < gamma_tmp)
		gamma0 = gamma_tmp;
	
	gamma.setConstant(gamma0+difference);
	
	x_bar_tmp = z.template segment<_n>((_pos_omega-1)*offset+_m+_n);	// depends on position of invariant set
	checkTheta = -f_xTheta + (F_xTheta*x_bar_tmp + F_theta*c_tmp);
	gamma0 = checkTheta.maxCoeff();
	gamma.template tail<1>().setConstant(gamma0+difference);
}


// ------------ function computes primal and dual residua needed by phaseI ---------------------
// stores in r_p (same as in PhaseII) and r_d_hat
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compRdRp_PhaseI()
{
	
	int offset1 = 2*(_n+1);	// offset required in nu for C'*nu
	
	// r_p = C*z - b_hat
	r_p_hat.template segment<_n>(0) = -Bm*z_hat.template segment<_m>(0) + z_hat.template segment<_n>(_m+1) - ( Am_tilde + Bm_bar)*x_hat - tm_tilde;
	r_p_hat.template segment<_n>(_n) = -B*z_hat.template segment<_m>(0) + z_hat.template segment<_n>(_m+1+_n) - (A_bar*x_hat + s);
	r_p_hat.template segment<1>(_n+_n) = z_hat.template segment<1>(_m) - z_hat.template segment<1>(_m+1+_n+_n);
	
	for (int i=1; i <= _N-1; i++)
	{
		r_p_hat.template segment<_n>((i-1)*offset1+_n+_n+1) = -Am_tilde*z_hat.template segment<_n>((i-1)*offset2+_m+1) - Bm_bar*z_hat.template segment<_n>((i-1)*offset2+_m+1+_n) -
			Bm*z_hat.template segment<_m>(i*offset2) + z_hat.template segment<_n>(i*offset2+_m+1) - tm_tilde;
		r_p_hat.template segment<1>((i-1)*2*(_n+1)+_n+_n+1+_n) = z_hat.template segment<1>((i-1)*offset2+_m+1+_n+_n) - 
			z_hat.template segment<1>(i*offset2+_m);

		r_p_hat.template segment<_n>(i*2*(_n+1)+_n) = -A_bar*z_hat.template segment<_n>((i-1)*offset2+_m+1+_n) - B*z_hat.template segment<_m>(i*offset2) + 
			z_hat.template segment<_n>(i*offset2+_m+1+_n) - s;
		r_p_hat.template segment<1>(i*2*(_n+1)+_n+_n) = z_hat.template segment<1>(i*offset2+_m) - z_hat.template segment<1>(i*offset2+_m+1+_n+_n);
	}
	// cout << setprecision(15) << "z_hat" << endl << z_hat << endl << endl;
	// cout << setprecision(15) << "r_p_hat in PhaseI:" << endl << r_p_hat << endl << endl;	
		
		
		
	// r_d_hat = 2*H_hat*z_hat + g_hat + kappa*P_hat'*d + C_hat'*nu_hat
	// note: H_hat = reg_hat * Identity(), g_hat has only two non-zero component
	Matrix<Type, _m+1, 1> g_hat_tmp;	// only for first and last block
	g_hat_tmp.setZero();
	g_hat_tmp.template tail<1>().setConstant(weight_PhaseI);
		
	r_d_hat.template segment<_m+1>(0) = 2*reg_hat*z_hat.template segment<_m+1>(0) + g_hat_tmp + 
			kappa * Fu_hat_transp[0] * d.template segment<_nInp>(0) - Bm_hat0_transp* nu_hat.template segment<_n>(0) - B_hat_transp * nu_hat.template segment<_n+1>(_n);		
		
	{	// do first step manually
		r_d_hat.template segment<_n>(_m+1) = 2*reg_hat*z_hat.template segment<_n>(_m+1) +
			nu_hat.template segment<_n>(0) - Am_tilde_transp*nu_hat.template segment<_n>(_n+_n+1);
		
		if (1 != _pos_omega)
		{
			r_d_hat.template segment<_n+1>(_m+1+_n) = 2*reg_hat*z_hat.template segment<_n+1>(_m+1+_n) +
				kappa*Fx_hat_transp[0] * d.template segment<_nSt>(_nInp) + kappa*Fu_bar_hat_transp[0]*d.template segment<_nInp>(_nInp+_nSt) + 
				Identity2_hat_transp*nu_hat.template segment<_n+1>(_n) - Bm_bar_hat_transp*nu_hat.template segment<_n+1>(_n+_n+1) - A_bar_hat_transp*nu_hat.template segment<_n+1>(_n+_n+1+_n+1);		
			// cout << setprecision(15) << "d in PhaseI" << endl << d << endl << endl;
			// cout << "nu_hat" << endl << nu_hat << endl << endl;
		}
		else	// 1 == _pos_omega
		{
			r_d_hat.template segment<_n+1>(_m+1+_n) = 2*reg_hat*z_hat.template segment<_n+1>(_m+1+_n) +
				kappa*Fx_hat_transp[0] * d.template segment<_nSt>(_nInp) + kappa*Fu_bar_hat_transp[0]*d.template segment<_nInp>(_nInp+_nSt) + 
				Identity2_hat_transp*nu_hat.template segment<_n+1>(_n) - Bm_bar_hat_transp*nu_hat.template segment<_n+1>(_n+_n+1) - A_bar_hat_transp*nu_hat.template segment<_n+1>(_n+_n+1+_n+1)
				+ kappa*F_xTheta_hat_transp*d.template segment<_nF_xTheta>(_N*(_nSt + _nInp));
		}
		r_d_hat.template segment<_m+1>(offset2) = 2*reg_hat*z_hat.template segment<_m+1>(offset2) +
			kappa*Fu_hat_transp[1]*d.template segment<_nInp>(_nInp + _nSt) - 
			Bm_hat_transp*nu_hat.template segment<_n+1>(_n+_n+1) - B_hat_transp*nu_hat.template segment<_n+1>(_n+_n+1+_n+1);		
	}
		
	for (int i=2; i<= _N-1; i++)
	{
		// note: H_hat = reg_hat * Identity();
		r_d_hat.template segment<_n>(_m+1+(i-1)*offset2) = 2*reg_hat*z_hat.template segment<_n>(_m+1+(i-1)*offset2) +
				nu_hat.template segment<_n>((i-2)*offset1+_n+_n+1) - Am_tilde_transp*nu_hat.template segment<_n>((i-1)*offset1+_n+_n+1);
			
		if (i != _pos_omega)
		{
			r_d_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = 2*reg_hat*z_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) +
					kappa*Fx_hat_transp[i-1] * d.template segment<_nSt>((i-1)*(_nInp+_nSt)+_nInp) + kappa*Fu_bar_hat_transp[i-1]*d.template segment<_nInp>(i*(_nSt+_nInp)) + 
					Identity2_hat_transp*nu_hat.template segment<_n+1>((i-1)*offset1+_n) - Bm_bar_hat_transp*nu_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) - A_bar_hat_transp*nu_hat.template segment<_n+1>(i*offset1+_n);		
		}
		else	// i == _pos_omega
		{
			r_d_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = 2*reg_hat*z_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) +
					kappa*Fx_hat_transp[i-1] * d.template segment<_nSt>((i-1)*(_nInp+_nSt)+_nInp) + kappa*Fu_bar_hat_transp[i-1]*d.template segment<_nInp>(i*(_nSt+_nInp)) + 
					Identity2_hat_transp*nu_hat.template segment<_n+1>((i-1)*offset1+_n) - Bm_bar_hat_transp*nu_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) - A_bar_hat_transp*nu_hat.template segment<_n+1>(i*offset1+_n) 
			 		+ kappa*F_xTheta_hat_transp*d.template segment<_nF_xTheta>(_N*(_nSt + _nInp));
		}
		r_d_hat.template segment<_m+1>(i*offset2) = 2*reg_hat*z_hat.template segment<_m+1>(i*offset2)  +
				kappa*Fu_hat_transp[i]*d.template segment<_nInp>(i*(_nInp+_nSt)) - 
				Bm_hat_transp*nu_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) - B_hat_transp*nu_hat.template segment<_n+1>(i*offset1+_n);		
	}
		
	r_d_hat.template segment<_n>(_m+1+(_N-1)*offset2) = 2*reg_hat*z_hat.template segment<_n>(_m+1+(_N-1)*offset2) + 
			 	nu_hat.template segment<_n>((_N-2)*offset1+_n+_n+1);
	if (_pos_omega == _N)
	{
		r_d_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = 2*reg_hat*z_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) + 
				kappa*Fx_hat_transp[_N-1]*d.template segment<_nSt>((_N-1)*(_nSt+_nInp)+_nInp) + kappa*F_xTheta_hat_transp*d.template segment<_nF_xTheta>(_N*(_nSt+_nInp)) +
				Identity2_hat_transp*nu_hat.template segment<_n+1>((_N-1)*offset1+_n);
	}
	else
	{
		r_d_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = 2*reg_hat*z_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) + 
				kappa*Fx_hat_transp[_N-1]*d.template segment<_nSt>((_N-1)*(_nSt+_nInp)+_nInp) + 
				Identity2_hat_transp*nu_hat.template segment<_n+1>((_N-1)*offset1+_n);
	}
	r_d_hat.template segment<_m+1>(_N*offset2) = 2*reg_hat*z_hat.template segment<_m+1>(_N*offset2) + g_hat_tmp + 
				kappa * F_theta_hat_transp*d.template segment<_nF_xTheta>(_N*(_nSt+_nInp));
	
	// cout << setprecision(30) << "r_d_hat:" << endl << r_d_hat << endl << endl;
}



// ---------------- computes the updated for z_hat and nu
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDz_hatDnu_hat_PhaseI()
{
	compPhi_hat_PhaseI();
	compPhi_hat_tilde_PhaseI();	// Phi_hat_tilde = Phi_hat^{-1}
	compY_hat_PhaseI();	// stored in Y_hat
	compBeta_hat_PhaseI();	// stored in beta_hat
	compL_hat();	// elements stored in L from PhaseII
	compDnu_hat();	// elements stored in dnu from PhaseII, Y*dnu = -beta_hat;
	compDz_hat_PhaseI();
}



// ------- computes Phi_hat = 2*H_hat + kappa*P'*diag(d)^2*P ------------------
// ------------ works ------------------------------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compPhi_hat_PhaseI()
{
	// ------------ first, partition diag(d) into (2N+1) diagonal matrices called d_diag[2N+1]
	
	/*
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
	*/
	
	DiagonalMatrix<Type, _nSt> d_diagX[_N];
	DiagonalMatrix<Type, _nInp> d_diagU[_N];
	DiagonalMatrix<Type, _nF_xTheta> d_diagTheta;
	int tmp = 0;
	for (int i=0 ; i <= _N-1 ; i++)
	{
		d_diagU[i].diagonal() = d.template segment<_nInp>(tmp);
		tmp = tmp + _nInp;
		d_diagX[i].diagonal() = d.template segment<_nSt>(tmp);
		tmp = tmp + _nSt;
	}
	d_diagTheta.diagonal() = d.template segment<_nF_xTheta>(tmp);
	
	
	Matrix<Type, _m+1, _m+1> eye_m1;
	eye_m1.setIdentity();
	Matrix<Type, _n+1, _n+1> eye_n1;
	eye_n1.setIdentity();

	// ------- compute elements of Phi_hat, where Pi_tilde from PhaseII can be used
	// Omicron_hat[0] = 2*reg_hat*eye_m1 + kappa*Fu_hat_transp[0]*d_diag[0]*d_diag[0]*Fu_hat[0];
	Omicron_hat[0] = 2*reg_hat*eye_m1 + kappa*Fu_hat_transp[0]*d_diagU[0]*d_diagU[0]*Fu_hat[0];
	
	for (int i=1; i <= _N-1; i++)
	{
		if (i != _pos_omega)
		{
			// Rho_hat[i-1] = 2*reg_hat*eye_n1  + kappa * ( Fx_hat_transp[i-1]*d_diag[2*(i-1)+1]*d_diag[2*(i-1)+1]*Fx_hat[i-1]  + Fu_bar_hat_transp[i-1]*d_diag[2*i]*d_diag[2*i]*Fu_bar_hat[i-1] );
			Rho_hat[i-1] = 2*reg_hat*eye_n1  + kappa * ( Fx_hat_transp[i-1]*d_diagX[i-1]*d_diagX[i-1]*Fx_hat[i-1]  + Fu_bar_hat_transp[i-1]*d_diagU[i]*d_diagU[i]*Fu_bar_hat[i-1] );
		}
		else	// i == _pos_omega
		{
			// Rho_hat[_pos_omega-1] = 2*reg_hat*eye_n1  + kappa * ( Fx_hat_transp[_pos_omega-1]*d_diag[2*(_pos_omega-1)+1]*d_diag[2*(_pos_omega-1)+1]*Fx_hat[_pos_omega-1]  + Fu_bar_hat_transp[_pos_omega-1]*d_diag[2*_pos_omega]*d_diag[2*_pos_omega]*Fu_bar_hat[_pos_omega-1] 
				// + F_xTheta_hat_transp * d_diag[2*_N]*d_diag[2*_N]*F_xTheta_hat );
			Rho_hat[_pos_omega-1] = 2*reg_hat*eye_n1  + kappa * ( Fx_hat_transp[_pos_omega-1]*d_diagX[_pos_omega-1]*d_diagX[_pos_omega-1]*Fx_hat[_pos_omega-1]  + Fu_bar_hat_transp[_pos_omega-1]*d_diagU[_pos_omega]*d_diagU[_pos_omega]*Fu_bar_hat[_pos_omega-1] 
				+ F_xTheta_hat_transp * d_diagTheta*d_diagTheta*F_xTheta_hat );
		}
		// Sigma_hat[i-1] = kappa * ( Fu_bar_hat_transp[i-1]*d_diag[2*i]*d_diag[2*i]*Fu_hat[i] );
		Sigma_hat[i-1] = kappa * ( Fu_bar_hat_transp[i-1]*d_diagU[i]*d_diagU[i]*Fu_hat[i] );
		// Omicron_hat[i] = 2*reg_hat*eye_m1 + kappa * ( Fu_hat_transp[i]*d_diag[2*i]*d_diag[2*i]*Fu_hat[i] );
		Omicron_hat[i] = 2*reg_hat*eye_m1 + kappa * ( Fu_hat_transp[i]*d_diagU[i]*d_diagU[i]*Fu_hat[i] );
	}
	
	// special treatment for last block
	if (_pos_omega == _N)
	{
		// Rho_hat[_N-1] = 2*reg_hat*eye_n1 + kappa * ( Fx_hat_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx_hat[_N-1] + F_xTheta_hat_transp*d_diag[2*_N]*d_diag[2*_N]*F_xTheta_hat );
		Rho_hat[_N-1] = 2*reg_hat*eye_n1 + kappa * ( Fx_hat_transp[_N-1]*d_diagX[_N-1]*d_diagX[_N-1]*Fx_hat[_N-1] + F_xTheta_hat_transp*d_diagTheta*d_diagTheta*F_xTheta_hat );
	}
	else
	{	
		// Rho_hat[_N-1] = 2*reg_hat*eye_n1 + kappa * ( Fx_hat_transp[_N-1]*d_diag[2*_N-1]*d_diag[2*_N-1]*Fx_hat[_N-1] );
		Rho_hat[_N-1] = 2*reg_hat*eye_n1 + kappa * ( Fx_hat_transp[_N-1]*d_diagX[_N-1]*d_diagX[_N-1]*Fx_hat[_N-1] );
	}
	
	// Sigma_hat[_N-1] = kappa * ( F_xTheta_hat_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta_hat );
	Sigma_hat[_N-1] = kappa * ( F_xTheta_hat_transp*d_diagTheta*d_diagTheta*F_theta_hat );
	// Omicron_hat[_N] = 2*reg_hat*eye_m1 + kappa * ( F_theta_hat_transp*d_diag[2*_N]*d_diag[2*_N]*F_theta_hat );
	Omicron_hat[_N] = 2*reg_hat*eye_m1 + kappa * ( F_theta_hat_transp*d_diagTheta*d_diagTheta*F_theta_hat );
	
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
	LOmicron_diag[0].compute(Omicron_hat[0]);
	LOmicron_hat_diag_transp[0] = LOmicron_diag[0].matrixLLT().transpose();
	Matrix<Type, _n, _n> eye;
	eye.setIdentity();
	
	for (int i=1 ; i<= _N-1; i++)
	{
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
			LLambda0_hat_transp = LRho_diag[_pos_omega-1].matrixLLT().triangularView<Lower>().solve(Sigma_hat[_N-1]);
			LLambda0_hat = LLambda0_hat_transp.transpose();
			LLambda1_hat_transp = LOmicron_diag[_pos_omega].matrixLLT().triangularView<Lower>().solve(Matrix<Type,_m+1,_m+1>::Zero() - LSigma_hat_offDiag[_pos_omega-1]*LLambda0_hat_transp);
			LLambda1_hat = LLambda1_hat_transp.transpose();
		}
	}
	
	LPi_diag[_N-1].compute(2*reg_hat*eye);
	LPi_hat_diag_transp[_N-1] = LPi_diag[_N-1].matrixLLT().transpose();
	LRho_diag[_N-1].compute(Rho_hat[_N-1]);
	LRho_hat_diag_transp[_N-1] = LRho_diag[_N-1].matrixLLT().transpose();		
	
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
	
	// 1. Compute elements of Matrix U_hat
	// treat the first 2 U_bar and first 3 U specially
	U_hat[0] = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(-Bm_hat0_transp);
	U_bar_hat[0] = LPi_diag[0].matrixLLT().triangularView<Lower>().solve(eye);
	
	U_hat[1] = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(-B_hat_transp);
	U_bar_hat[1] = LRho_diag[0].matrixLLT().triangularView<Lower>().solve(Identity2_hat_transp);
	U_hat[2] = LOmicron_diag[1].matrixLLT().triangularView<Lower>().solve( - LSigma_hat_offDiag[0]*U_bar_hat[1] );
	
	// remaining U_bar and U have structures
	for (int i=1; i<= _N-2; i++)
	{
		U_bar_hat[2+(i-1)*5] = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(-Am_tilde_hat_transp);
		U_bar_hat[3+(i-1)*5] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(-Bm_bar_hat_transp);
		U_hat[3+(i-1)*3] = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(-Bm_hat_transp - LSigma_hat_offDiag[i-1]*U_bar_hat[3+(i-1)*5] );
		U_bar_hat[4+(i-1)*5] = LPi_diag[i].matrixLLT().triangularView<Lower>().solve( Identity1_hat_transp );
		
		U_bar_hat[5+(i-1)*5] = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(-A_bar_hat_transp);
		U_hat[4+(i-1)*3] = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(-B_hat_transp - LSigma_hat_offDiag[i-1]*U_bar_hat[5+(i-1)*5]);
		U_bar_hat[6+(i-1)*5] = LRho_diag[i].matrixLLT().triangularView<Lower>().solve( Identity2_hat_transp );
		U_hat[5+(i-1)*3] = LOmicron_diag[i+1].matrixLLT().triangularView<Lower>().solve( -LSigma_hat_offDiag[i]*U_bar_hat[6+(i-1)*5] );
	}
	
	
	U_bar_hat[2+(_N-2)*5] = LPi_diag[_N-2].matrixLLT().triangularView<Lower>().solve(-Am_tilde_hat_transp);
	U_bar_hat[3+(_N-2)*5] = LRho_diag[_N-2].matrixLLT().triangularView<Lower>().solve(-Bm_bar_hat_transp);
	U_hat[3+(_N-2)*3] = LOmicron_diag[_N-1].matrixLLT().triangularView<Lower>().solve(-Bm_hat_transp - LSigma_hat_offDiag[_N-2]*U_bar_hat[3+(_N-2)*5] );
	U_bar_hat[4+(_N-2)*5] = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve( Identity1_hat_transp );
	
	U_bar_hat[5+(_N-2)*5] = LRho_diag[_N-1-1].matrixLLT().triangularView<Lower>().solve(-A_bar_hat_transp);
	U_hat[4+(_N-2)*3] = LOmicron_diag[_N-1].matrixLLT().triangularView<Lower>().solve(-B_hat_transp - LSigma_hat_offDiag[_N-1-1]*U_bar_hat[5+(_N-1-1)*5]);
	U_bar_hat[6+(_N-2)*5] = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve( Identity2_hat_transp );
	
	if (_N == _pos_omega)
	{
		U_hat[5+(_N-2)*3] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( -LSigma_hat_offDiag[_N-1]*U_bar_hat[6+(_N-1-1)*5] );
	}
	else
	{
		UO_hat[0] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( -LLambda0_hat*U_bar_hat[1+(_pos_omega-1)*5] - LLambda1_hat*U_hat[2+(_pos_omega-1)*3]);
		UO_hat[1] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( -LLambda0_hat*U_bar_hat[3+(_pos_omega-1)*5] - LLambda1_hat*U_hat[3+(_pos_omega-1)*3]);
		UO_hat[2] = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve( -LLambda0_hat*U_bar_hat[5+(_pos_omega-1)*5] - LLambda1_hat*U_hat[4+(_pos_omega-1)*3]);
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
	for (int i = 1; i <= _N-2; i++)
	{
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
	Y_hat[0][0] = -Bm_hat0*X_hat[0] + X_bar_hat[0];
	
	Y_hat[0][1] = -Bm_hat0*X_hat[1];
	Y_hat[1][1] = -B_hat*X_hat[1] + Identity2_hat*X_bar_hat[1];
	
	{	// do another one manually b/c first block-row of Y has a row less
		Y_hat[0][2] = X_bar_hat[2];
		Y_hat[1][2] = Identity2_hat*X_bar_hat[3];
		Y_hat[2][2] = -Am_tilde_hat*X_bar_hat[2] - Bm_bar_hat*X_bar_hat[3] - Bm_hat*X_hat[3] + Identity1_hat*X_bar_hat[4];
		
		Y_hat[0][3] = Identity2_hat*X_bar_hat[5];
		Y_hat[1][3] = -Bm_bar_hat*X_bar_hat[5] - Bm_hat*X_hat[4];
		Y_hat[2][3] = -A_bar_hat*X_bar_hat[5] - B_hat*X_hat[4] + Identity2_hat*X_bar_hat[6];
	}
	
	
	// compute rest by filling Y column by column; done for 2 neighboring rows
	// we fill Y column by column, treating the first two columns specially
	for (int i=2; i <= _N-1; i++)
	{
		Y_hat[0][2*(i-1)+2] = Identity1_hat*X_bar_hat[2+(i-1)*5];
		Y_hat[1][2*(i-1)+2] = Identity2_hat*X_bar_hat[3+(i-1)*5];
		Y_hat[2][2*(i-1)+2] = -Am_tilde_hat*X_bar_hat[2+(i-1)*5] - Bm_bar_hat*X_bar_hat[3+(i-1)*5] - Bm_hat*X_hat[3+(i-1)*3] + Identity1_hat*X_bar_hat[4+(i-1)*5];
			
		Y_hat[0][2*(i-1)+3] = Identity2_hat*X_bar_hat[5+(i-1)*5];
		Y_hat[1][2*(i-1)+3] = -Bm_bar_hat*X_bar_hat[5+(i-1)*5] - Bm_hat*X_hat[4+(i-1)*3];
		Y_hat[2][2*(i-1)+3] = -A_bar_hat*X_bar_hat[5+(i-1)*5] - B_hat*X_hat[4+(i-1)*3] + Identity2_hat*X_bar_hat[6+(i-1)*5];
	}
	
	/*
	cout << "Y_hat[0][0]" << endl << Y_hat[0][0] << endl << endl; 
	cout << "Y_hat[0][1]" << endl << Y_hat[0][1] << endl << endl;
	cout << "Y_hat[1][1]" << endl << Y_hat[1][1] << endl << endl;;
			
	for (int i=1; i <= _N-1; i++)
	{
		cout << "Y_hat[0][" << 2*(i-1)+3 << "]" << endl << Y_hat[0][2*(i-1)+2] << endl << endl; 
		cout << " Y_hat[1][" << 2*(i-1)+3 << "]" << endl << Y_hat[1][2*(i-1)+2] << endl << endl;
		cout << "Y_hat[2][" << 2*(i-1)+3 << "]" << endl << Y_hat[2][2*(i-1)+2] << endl << endl;
				
		cout << "Y_hat[0][" << 2*(i-1)+3 << "]" << endl << Y_hat[0][2*(i-1)+3] << endl << endl; 
		cout << "Y_hat[1][" << 2*(i-1)+3 << "]" << endl << Y_hat[1][2*(i-1)+3] << endl << endl; 
		cout << "Y_hat[2][" << 2*(i-1)+3 << "]" << endl << Y_hat[2][2*(i-1)+3] << endl << endl; 
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
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> tmp2_hat;	// long vector
	tmp2_hat.template segment<_m+1>(0) = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+1>(0) );
	
	for (int i = 1; i <= _N-1; i++)
	{
		tmp2_hat.template segment<_n>(_m+1+(i-1)*offset2) = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n>(_m+1+(i-1)*offset2));
		tmp2_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n));
		tmp2_hat.template segment<_m+1>(i*offset2) = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+1>(i*offset2) - LSigma_hat_offDiag[i-1]*tmp2_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n));
	}
	
	
	if(_pos_omega == _N)
	{	
		tmp2_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n>(_m+1+(_N-1)*offset2));
		tmp2_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
		tmp2_hat.template segment<_m+1>(_N*offset2) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+1>(_N*offset2) - LSigma_hat_offDiag[_N-1]*tmp2_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
	}
	else	// use LLambda0_hat and LLambda1_hat
	{
		tmp2_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n>(_m+1+(_N-1)*offset2));
		tmp2_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
		tmp2_hat.template segment<_m+1>(_N*offset2) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(r_d_hat.template segment<_m+1>(_N*offset2) - LLambda0_hat*tmp2_hat.template segment<_n+1>(_m+1+(_pos_omega-1)*offset2+_n) - LLambda1_hat*tmp2_hat.template segment<_m+1>(_m+1+(_pos_omega-1)*offset2+_n+_n+1)   );
	}
	// cout << setprecision(20) << "tmp2_hat:" << endl << tmp2_hat << endl << endl;
	
	// 2. compute tmp1: L'*tmp1 = tmp2
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> tmp1_hat;	// long vector
	tmp1_hat.template segment<_m+1>(0) = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+1>(0));
	
	for (int i = 1; i <= _pos_omega-1; i++)
	{
		tmp1_hat.template segment<_n>(_m+1+(i-1)*offset2) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+1+(i-1)*offset2));
		tmp1_hat.template segment<_m+1>(i*offset2) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+1>(i*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) - LSigma_hat_offDiag_transp[i-1]*tmp1_hat.template segment<_m+1>(i*offset2));
	}
	
	// the missing block is computed after the last block is computed
	for (int i = _pos_omega+1; i <= _N-1; i++)
	{
		tmp1_hat.template segment<_n>(_m+1+(i-1)*offset2) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+1+(i-1)*offset2));
		tmp1_hat.template segment<_m+1>(i*offset2) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+1>(i*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) - LSigma_hat_offDiag_transp[i-1]*tmp1_hat.template segment<_m+1>(i*offset2));
	}
	
	if (_pos_omega == _N)
	{
		tmp1_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+1+(_pos_omega-1)*offset2));
		tmp1_hat.template segment<_m+1>(_N*offset2) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+1>(_N*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) - LSigma_hat_offDiag_transp[_N-1]*tmp1_hat.template segment<_m+1>(_N*offset2));
	}
	else	// standard ending, compute the missing blocks
	{
		tmp1_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+1+(_N-1)*offset2));
		tmp1_hat.template segment<_m+1>(_N*offset2) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+1>(_N*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
		
		tmp1_hat.template segment<_n>(_m+1+(_pos_omega-1)*offset2) = LPi_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n>(_m+1+(_pos_omega-1)*offset2) );
		tmp1_hat.template segment<_m+1>(_pos_omega*offset2) = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp2_hat.template segment<_m+1>(_pos_omega*offset2) - LLambda1_hat_transp*tmp1_hat.template segment<_m+1>(_N*offset2)  );
		tmp1_hat.template segment<_n+1>(_m+1+(_pos_omega-1)*offset2+_n) = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp2_hat.template segment<_n+1>(_m+1+(_pos_omega-1)*offset2+_n) - LSigma_hat_offDiag_transp[_pos_omega-1]*tmp1_hat.template segment<_m+1>(_pos_omega*offset2) - LLambda0_hat_transp*tmp1_hat.template segment<_m+1>(_N*offset2)  );
	}
	// cout << setprecision(30) << "tmp1_hat:" << endl << tmp1_hat << endl << endl;
	
	// 3. beta_hat = -r_p_hat + C_hat * tmp1_hat
	beta_hat.template segment<_n>(0) = -r_p_hat.template segment<_n>(0) + ( -Bm*tmp1_hat.template segment<_m>(0) + tmp1_hat.template segment<_n>(_m+1) );
	beta_hat.template segment<_n+1>(_n) = -r_p_hat.template segment<_n+1>(_n) + (- B_hat*tmp1_hat.template segment<_m+1>(0) + Identity2_hat*tmp1_hat.template segment<_n+1>(_m+1+_n) );
	
	int offset1 = 2*(_n+1);
	for(int i=1; i<= _N-1; i++)
	{
		beta_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) = -r_p_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) + (    
			- Am_tilde_hat * tmp1_hat.template segment<_n>(_m+1+(i-1)*offset2) 
			- Bm_bar_hat * tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n)
			- Bm_hat * tmp1_hat.template segment<_m+1>(i*offset2) 
			+ Identity1_hat * tmp1_hat.template segment<_n>(i*offset2+_m+1)    );
					
		beta_hat.template segment<_n+1>(i*offset1+_n) = -r_p_hat.template segment<_n+1>( i*offset1+_n ) + (  
		    - A_bar_hat * tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n)
		    - B_hat * tmp1_hat.template segment<_m+1>(i*offset2)
			+ Identity2_hat * tmp1_hat.template segment<_n+1>(i*offset2+_m+1+_n )        );
	}
	// cout << setprecision(15) << "beta_hat PhaseI:" << endl << beta_hat << endl << endl;	
}

// ------------ function computes L_hat: L_hat*L_hat' = Y_hat --------------------------
// -------------- compute components of L[3][2*_N] ---------------------
// ---- remark: L[i][i]=L_diag are lower triangular matrices and can be accessed by L[i][i].matrixLLT().triangularView<Lower>()
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compL_hat()		// L*L' = Y
{
	// special treatment for L11, L21, L31
	L_diag[0].compute(Y_hat[0][0]);		// Y11 = L11*L11'
	L_hat_diag_transp[0] = L_diag[0].matrixLLT().transpose();
	
	L_hat_offDiag_transp[0][0] = L_diag[0].matrixLLT().triangularView<Lower>().solve(Y_hat[0][1]);
	L_hat_offDiag[0][0] = L_hat_offDiag_transp[0][0].transpose();
	
	L_hat_offDiag_transp[1][0] = L_diag[0].matrixLLT().triangularView<Lower>().solve(Y_hat[0][2]);
	L_hat_offDiag[1][0] = L_hat_offDiag_transp[1][0].transpose();	
	
	
	// special treatment for L22, L32, L42
	L_diag[1].compute( Y_hat[1][1]-L_hat_offDiag[0][0]*L_hat_offDiag_transp[0][0] );
	L_hat_diag_transp[1] = L_diag[1].matrixLLT().transpose();
	
	L_hat_offDiag_transp[0][1] = L_diag[1].matrixLLT().triangularView<Lower>().solve(Y_hat[1][2]-L_hat_offDiag[0][0]*L_hat_offDiag_transp[1][0]);
	L_hat_offDiag[0][1] = L_hat_offDiag_transp[0][1].transpose();
	
	L_hat_offDiag_transp[1][1] = L_diag[1].matrixLLT().triangularView<Lower>().solve(Y_hat[0][3]);
	L_hat_offDiag[1][1] = L_hat_offDiag_transp[1][1].transpose();
	
	
	// cases in the middle
	for (int i = 1; i <= 2*_N-4; i++)
	{
		L_diag[i+1].compute( Y_hat[2][i+1] - L_hat_offDiag[1][i-1]*L_hat_offDiag_transp[1][i-1] - L_hat_offDiag[0][i]*L_hat_offDiag_transp[0][i] );
		L_hat_diag_transp[i+1] = L_diag[i+1].matrixLLT().transpose();
		L_hat_offDiag_transp[0][i+1] = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( Y_hat[1][i+2] - L_hat_offDiag[0][i]*L_hat_offDiag_transp[1][i] );
		L_hat_offDiag[0][i+1] = L_hat_offDiag_transp[0][i+1].transpose();
		
		L_hat_offDiag_transp[1][i+1] = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( Y_hat[0][i+3] );
		L_hat_offDiag[1][i+1] = L_hat_offDiag_transp[1][i+1].transpose();
	}
			
	// special treatment in the end, i.e. i = 2*_N-3
	L_diag[2*_N-2].compute( Y_hat[2][2*_N-2] - L_hat_offDiag[1][2*_N-4]*L_hat_offDiag_transp[1][2*_N-4] - L_hat_offDiag[0][2*_N-3]*L_hat_offDiag_transp[0][2*_N-3] );
	L_hat_diag_transp[2*_N-2] = L_diag[2*_N-2].matrixLLT().transpose();
	L_hat_offDiag_transp[0][2*_N-2] = L_diag[2*_N-2].matrixLLT().triangularView<Lower>().solve( Y_hat[1][2*_N-1] - L_hat_offDiag[0][2*_N-3]*L_hat_offDiag_transp[1][2*_N-3] );
	L_hat_offDiag[0][2*_N-2] = L_hat_offDiag_transp[0][2*_N-2].transpose();
	
	// i = 2*_N-2
	L_diag[2*_N-1].compute( Y_hat[2][2*_N-1] - L_hat_offDiag[1][2*_N-3]*L_hat_offDiag_transp[1][2*_N-3] - L_hat_offDiag[0][2*_N-2]*L_hat_offDiag_transp[0][2*_N-2] );
	L_hat_diag_transp[2*_N-1] = L_diag[2*_N-1].matrixLLT().transpose();	
}


// ------------ function computes L_hat*L_hat'*dnu_hat = -beta_hat ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDnu_hat()
{
	
	// 1) first, solve for delta: L_hat*delta_hat = -beta_hat
	Matrix<Type, 2*_N*(_n+1)-1, 1> delta_hat;
	
	// special cases in the beginning
	delta_hat.template segment<_n>(0) = L_diag[0].matrixLLT().triangularView<Lower>().solve(-beta_hat.template segment<_n>(0));
	delta_hat.template segment<_n+1>(_n) = L_diag[1].matrixLLT().triangularView<Lower>().solve(-beta_hat.template segment<_n+1>(_n) - L_hat_offDiag[0][0]*delta_hat.template segment<_n>(0));
	delta_hat.template segment<_n+1>(_n+_n+1) = L_diag[2].matrixLLT().triangularView<Lower>().solve( -beta_hat.template segment<_n+1>(_n+_n+1) - L_hat_offDiag[1][0]*delta_hat.template segment<_n>(0) - L_hat_offDiag[0][1]*delta_hat.template segment<_n+1>(_n) );
	
	// remaining cases are regular
	for (int i=2; i<= 2*_N-2; i++)
	{
		delta_hat.template segment<_n+1>(_n+i*(_n+1)) = L_diag[i+1].matrixLLT().triangularView<Lower>().solve( -beta_hat.template segment<_n+1>(_n+i*(_n+1)) - L_hat_offDiag[1][i-1]*delta_hat.template segment<_n+1>((i-2)*(_n+1)+_n) - L_hat_offDiag[0][i]*delta_hat.template segment<_n+1>((i-1)*(_n+1)+_n) );
	}
	// cout << setprecision(30) << "delta_hat:" << endl << delta_hat << endl << endl;
	
	
	// 2) now, solve for L'*Dnu_hat = delta_hat
	// from behind...
	int length = 2*_N*(_n+1)-1;
	int tmp_offset = _n+1;
	dnu_hat.template segment<_n+1>(length - tmp_offset) = L_hat_diag_transp[2*_N-1].template triangularView<Upper>().solve( delta_hat.template segment<_n+1>(length - tmp_offset) );
	dnu_hat.template segment<_n+1>(length - 2*tmp_offset) = L_hat_diag_transp[2*_N-2].template triangularView<Upper>().solve( delta_hat.template segment<_n+1>(length - 2*tmp_offset) - L_hat_offDiag_transp[0][2*_N-2]*dnu_hat.template segment<_n+1>(length - tmp_offset) );
	
	//remaining cases are regular
	for (int i=1; i<=2*_N-3; i++)
	{
		dnu_hat.template segment<_n+1>(length-(i+2)*tmp_offset) = L_hat_diag_transp[2*_N-(i+2)].template triangularView<Upper>().solve( delta_hat.template segment<_n+1>(length-(i+2)*tmp_offset) - L_hat_offDiag_transp[0][2*_N-(i+2)]*dnu_hat.template segment<_n+1>(length-(i+1)*tmp_offset) - L_hat_offDiag_transp[1][2*_N-(i+2)]*dnu_hat.template segment<_n+1>(length-i*tmp_offset)  );
	}
	dnu_hat.template segment<_n>(0) = L_hat_diag_transp[0].template triangularView<Upper>().solve( delta_hat.template segment<_n>(0) - L_hat_offDiag_transp[0][0]*dnu_hat.template segment<_n+1>(_n) - L_hat_offDiag_transp[1][0]*dnu_hat.template segment<_n+1>(_n+_n+1)  );
	// cout << setprecision(30) << "dnu_hat" << endl << dnu_hat << endl << endl;
}



// ---- function computes Phi_hat * dz_hat = -r_d_hat - C_hat' * dnu_hat ----------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compDz_hat_PhaseI()
{
	
	// computed in parts
	// 1. tmp_hat = -r_d_hat - C_hat' * dnu
	// 2. L_hat*L_hat'*dz_hat = tmp_hat
	// 3. L_hat*tmp1_hat = tmp_hat
	// 4. L_hat'*dz_hat = tmp1_hat
	
	// 1. compute tmp_hat
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> tmp_hat;	// long vector
	
	tmp_hat.template segment<_m+1>(0) = -r_d_hat.template segment<_m+1>(0) + Bm_hat0_transp*dnu_hat.template segment<_n>(0) + B_hat_transp*dnu_hat.template segment<_n+1>(_n);
	int offset1 = 2*(_n+1);
	{	// do second block manually
		tmp_hat.template segment<_n>(_m+1) = -r_d_hat.template segment<_n>(_m+1) - dnu_hat.template segment<_n>(0) + Am_tilde_transp*dnu_hat.template segment<_n>(_n+_n+1) ;
		tmp_hat.template segment<_n+1>(_m+1+_n) = -r_d_hat.template segment<_n+1>(_m+1+_n) - Identity2_hat_transp*dnu_hat.template segment<_n+1>(_n) + Bm_bar_hat_transp*dnu_hat.template segment<_n+1>(_n+_n+1) + A_bar_hat_transp * dnu_hat.template segment<_n+1>(offset1+_n);
		tmp_hat.template segment<_m+1>(offset2) = -r_d_hat.template segment<_m+1>(offset2) + Bm_hat_transp*dnu_hat.template segment<_n+1>(_n+1+_n) + B_hat_transp*dnu_hat.template segment<_n+1>(offset1+_n);
	}
	// cout << setprecision(30) << "tmp_hat:" << endl << tmp_hat << endl << endl;
	
	for (int i=2; i<= _N-1; i++)
	{
		tmp_hat.template segment<_n>(_m+1+(i-1)*offset2) = -r_d_hat.template segment<_n>(_m+1+(i-1)*offset2) - dnu_hat.template segment<_n>((i-2)*offset1+_n+_n+1) + Am_tilde_transp*dnu_hat.template segment<_n>((i-1)*offset1+_n+_n+1) ;
		tmp_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = -r_d_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) - Identity2_hat_transp*dnu_hat.template segment<_n+1>((i-1)*offset1+_n) + Bm_bar_hat_transp*dnu_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) + A_bar_hat_transp * dnu_hat.template segment<_n+1>(i*offset1+_n);
		tmp_hat.template segment<_m+1>(i*offset2) = -r_d_hat.template segment<_m+1>(i*offset2) + Bm_hat_transp*dnu_hat.template segment<_n+1>((i-1)*offset1+_n+_n+1) + B_hat_transp*dnu_hat.template segment<_n+1>(i*offset1+_n);
	}
	
	
	tmp_hat.template segment<_n>(_m+1+(_N-1)*offset2) = -r_d_hat.template segment<_n>(_m+1+(_N-1)*offset2) - dnu_hat.template segment<_n>((_N-2)*offset1+_n+_n+1);
	tmp_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = -r_d_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) - Identity2_hat_transp*dnu_hat.template segment<_n+1>((_N-1)*offset1+_n);
	tmp_hat.template segment<_m+1>(_N*offset2) = -r_d_hat.template segment<_m+1>(_N*offset2);
	// cout << setprecision(30) << "tmp_hat:" << endl << tmp_hat << endl << endl;


	// 3. L_Phi*tmp1_hat = tmp_hat
	Matrix<Type, _N*(_m + _n + _n + 1 + 1) + _m + 1 , 1> tmp1_hat;	// long vector
	tmp1_hat.template segment<_m+1>(0) = LOmicron_diag[0].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+1>(0) );
	for (int i = 1; i <= _N-1; i++)
	{
		tmp1_hat.template segment<_n>(_m+1+(i-1)*offset2) = LPi_diag[i-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n>(_m+1+(i-1)*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = LRho_diag[i-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n));
		tmp1_hat.template segment<_m+1>(i*offset2) = LOmicron_diag[i].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+1>(i*offset2) - LSigma_hat_offDiag[i-1]*tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n));
	}
	
	if(_pos_omega == _N)
	{
		tmp1_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n>(_m+1+(_N-1)*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
		tmp1_hat.template segment<_m+1>(_N*offset2) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+1>(_N*offset2) - LSigma_hat_offDiag[_N-1]*tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
	}
	else
	{
		tmp1_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n>(_m+1+(_N-1)*offset2));
		tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_diag[_N-1].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));
		tmp1_hat.template segment<_m+1>(_N*offset2) = LOmicron_diag[_N].matrixLLT().triangularView<Lower>().solve(tmp_hat.template segment<_m+1>(_N*offset2) - LLambda0_hat*tmp1_hat.template segment<_n+1>(_m+1+(_pos_omega-1)*offset2+_n)  - LLambda1_hat*tmp1_hat.template segment<_m+1>(_m+1+(_pos_omega-1)*offset2+_n+_n+1) );
	}
	
	// cout << setprecision(30) << "tmp1_hat:" << endl << tmp1_hat << endl << endl;
	
	
	// 4. L_Phi'*dz_hat = tmp1_hat
	dz_hat.template segment<_m+1>(0) = LOmicron_hat_diag_transp[0].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+1>(0));
	for (int i = 1; i <= _pos_omega-1; i++)
	{
		dz_hat.template segment<_n>(_m+1+(i-1)*offset2) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+1+(i-1)*offset2));
		dz_hat.template segment<_m+1>(i*offset2) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+1>(i*offset2));
		dz_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) - LSigma_hat_offDiag_transp[i-1]*dz_hat.template segment<_m+1>(i*offset2) );
	}
	
	
	for (int i = _pos_omega+1; i <= _N-1; i++)
	{
		dz_hat.template segment<_n>(_m+1+(i-1)*offset2) = LPi_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+1+(i-1)*offset2));
		dz_hat.template segment<_m+1>(i*offset2) = LOmicron_hat_diag_transp[i].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+1>(i*offset2));
		dz_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) = LRho_hat_diag_transp[i-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+1>(_m+1+(i-1)*offset2+_n) - LSigma_hat_offDiag_transp[i-1]*dz_hat.template segment<_m+1>(i*offset2) );
	}
	
	
	if (_pos_omega == _N)
	{
		dz_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+1+(_pos_omega-1)*offset2));
		dz_hat.template segment<_m+1>(_N*offset2) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+1>(_N*offset2));
		dz_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) - LSigma_hat_offDiag_transp[_N-1]*dz_hat.template segment<_m+1>(_N*offset2) );
	}
	else 	// standard ending, compute missing block in middle
	{
		dz_hat.template segment<_n>(_m+1+(_N-1)*offset2) = LPi_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+1+(_N-1)*offset2));
		dz_hat.template segment<_m+1>(_N*offset2) = LOmicron_hat_diag_transp[_N].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+1>(_N*offset2));
		dz_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n) = LRho_hat_diag_transp[_N-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+1>(_m+1+(_N-1)*offset2+_n));

		dz_hat.template segment<_n>(_m+1+(_pos_omega-1)*offset2) = LPi_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n>(_m+1+(_pos_omega-1)*offset2) );
		dz_hat.template segment<_m+1>(_pos_omega*offset2) = LOmicron_hat_diag_transp[_pos_omega].template triangularView<Upper>().solve(tmp1_hat.template segment<_m+1>(_pos_omega*offset2) - LLambda1_hat_transp*dz_hat.template segment<_m+1>(_N*offset2)  );
		dz_hat.template segment<_n+1>(_m+1+(_pos_omega-1)*offset2+_n) = LRho_hat_diag_transp[_pos_omega-1].template triangularView<Upper>().solve(tmp1_hat.template segment<_n+1>(_m+1+(_pos_omega-1)*offset2+_n) - LSigma_hat_offDiag_transp[_pos_omega-1]*dz_hat.template segment<_m+1>(_pos_omega*offset2) - LLambda0_hat_transp*dz_hat.template segment<_m+1>(_N*offset2)  );
	}
	// cout << setprecision(30) << "dz_hat:" << endl << dz_hat << endl << endl;
}

/*
// ------------ function computes a feasibe gamma needed by phaseI ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: testPos_PhaseI()
{
	// special treatment at beginning
	c_hat_tmp = z_hat.template segment<_m+1>(0);
	checkU = (fu[0] - Fu[0]*K*x_hat) - (Fu_hat[0]*c_hat_tmp);	// should be >0
	for (int j=1; j <= _nInp; j++)
	{
		if (checkU[j-1] <= 0)
			return 0;
	}
		
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		x_bar_hat_tmp = z_hat.template segment<_n+1>((i-1)*offset2+_m+_n+1);
		c_hat_tmp = z_hat.template segment<_m+1>(i*offset2);		
		checkX = fx[i-1] - Fx_hat[i-1]*x_bar_hat_tmp;
		for (int j=1; j<= _nSt; j++)
		{
			if (checkX[j-1] <= 0)
				return 0;
		}
		checkU = fu[i] - (Fu_bar_hat[i-1]*x_bar_hat_tmp + Fu_hat[i]*c_hat_tmp);
		for (int j=1; j<=_nInp; j++)
		{
			if (checkU[j-1] <= 0)
				return 0;
		}
	}
	// cout << "finished general treatment at the beginning" << endl << endl;
		
	x_bar_hat_tmp = z_hat.template segment<_n+1>((_N-1)*offset2+_m+_n+1);
	c_hat_tmp = z_hat.template segment<_m+1>(_N*offset2);
	checkX = fx[_N-1] - (Fx_hat[_N-1]*x_bar_hat_tmp);
	for (int j=1; j<= _nSt; j++)
	{
		if (checkX[j-1] <= 0)
			return 0;
	}
	x_bar_hat_tmp = z_hat.template segment<_n+1>((_pos_omega-1)*offset2+_m+_n+1);
	checkTheta = f_xTheta - (F_xTheta_hat*x_bar_hat_tmp + F_theta_hat*c_hat_tmp);
	// cout << "check:" << endl << check << endl << endl;
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if (checkTheta[j-1] <= 0)
			return 0;
	}
	return 1;
}
*/


// ------------ function computes (d_hat)_i = (h - P_hat*z_hat)_i  --------
// elements are stored in d
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
void LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compD_hat_PhaseI()
{
	int tmp = 0;		// variable used to count position and offset
	
	// special treatment at beginning
	c_hat_tmp = z_hat.template segment<_m+1>(0);
	checkU = (fu[0] - Fu[0]*K*x_hat) - (Fu_hat[0]*c_hat_tmp);
	
	for (int j=1; j <= _nInp; j++)
		d[j-1] = 1/checkU[j-1];
	
	tmp = tmp + _nInp;
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		x_bar_hat_tmp = z_hat.template segment<_n+1>((i-1)*offset2+_m+_n+1);
		c_hat_tmp = z_hat.template segment<_m+1>(i*offset2);		
		checkX = fx[i-1] - Fx_hat[i-1]*x_bar_hat_tmp;
		for (int j=1; j<= _nSt; j++)
			d[tmp+j-1] = 1/checkX[j-1];
		
		tmp = tmp + _nSt;
		
		checkU = fu[i] - (Fu_bar_hat[i-1]*x_bar_hat_tmp + Fu_hat[i]*c_hat_tmp);
		for (int j=1; j<=_nInp; j++)
			d[tmp+j-1] = 1/checkU[j-1];
		
		tmp = tmp + _nInp;
	}
	
	x_bar_hat_tmp = z_hat.template segment<_n+1>((_N-1)*offset2+_m+_n+1);
	c_hat_tmp = z_hat.template segment<_m+1>(_N*offset2);
	checkX = fx[_N-1] - (Fx_hat[_N-1]*x_bar_hat_tmp);
	for (int j=1; j<= _nSt; j++)
		d[tmp+j-1] = 1/checkX[j-1];
	
	tmp = tmp + _nSt;
	x_bar_hat_tmp = z_hat.template segment<_n+1>((_pos_omega-1)*offset2+_m+_n+1);
	checkTheta = f_xTheta - (F_xTheta_hat*x_bar_hat_tmp + F_theta_hat*c_hat_tmp);
	for (int j=1; j<= _nF_xTheta; j++)
		d[tmp+j-1] = 1/checkTheta[j-1];
	
	// cout << setprecision(15) << "d in PhaseI" << endl << d << endl << endl;
	return;
}




// ------------ function checks if gamma < 0 and extracts z_warm  --------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
bool LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: testNeg_PhaseI()
{
	
	Matrix<Type, 2, 1> check;		// stores gamma segments
	
	for (int i=1; i <= _N; i++) // blocks in the middle
	{	
		check.template segment<1>(0) = z_hat.template segment<1>(_m+(i-1)*offset2);
		check.template segment<1>(1) = z_hat.template segment<1>(_m+1+_n+_n+(i-1)*offset2);
		// cout << "check" << endl << check << endl << endl;
		if( check[0] >= 0 || check[1] >= 0  )
			return 0;
		else
		{
			z_warm.template segment<_m>((i-1)*offset) = z_hat.template segment<_m>((i-1)*offset2);
			z_warm.template segment<_n+_n>((i-1)*offset + _m) = z_hat.template segment<_n+_n>((i-1)*offset2+_m+1);
		}
	}
	
	
	check.template segment<1>(0) = z_hat.template segment<1>(_N*offset2 + _m);
	if (check[0] >= 0)
		return 0;
	else
		z_warm.template segment<_m>(_N*offset ) = z_hat.template segment<_m>(_N*offset2);
	return 1;
}




// ------------ function returns largest t: P_hat*(z_hat+t*dz_hat) < h ---------------------
template <class Type, int _n, int _m, int _N, int _nSt, int _nInp, int _nF_xTheta, int _pos_omega>
double LBmpcTP<Type, _n, _m, _N, _nSt, _nInp, _nF_xTheta, _pos_omega> :: compT_PhaseI()
{
	
	double t=1;		// stores the smallest variable
	double t_tmp;	// t_tmp = min(t_vec)
	
	c_hat_tmp = z_hat.template segment<_m+1>(0);
	dc_hat_tmp = dz_hat.template segment<_m+1>(0);
	checkU = (fu[0] - Fu[0]*K*x_hat) - (Fu_hat[0]*c_hat_tmp);	// should be >0
	dcheckU = Fu_hat[0]*dc_hat_tmp;
	t_vecU.setConstant(1);
	for (int j=1; j <= _nInp; j++)
	{
		if (dcheckU[j-1] > 0)	// neg. cases not interesting
			t_vecU[j-1] = checkU[j-1]/dcheckU[j-1];
	}
	t_tmp = t_vecU.minCoeff();
	if (t_tmp < t)
		t = t_tmp;
	
	// general treatment in the middle; offset = _n + _n + _m;
	for (int i=1; i <= _N-1; i++) // blocks in the middle
	{	
		x_bar_hat_tmp = z_hat.template segment<_n+1>((i-1)*offset2+_m+_n+1);
		dx_bar_hat_tmp = dz_hat.template segment<_n+1>((i-1)*offset2+_m+_n+1);
		c_hat_tmp = z_hat.template segment<_m+1>(i*offset2);
		dc_hat_tmp = dz_hat.template segment<_m+1>(i*offset2);		
		checkX = fx[i-1] - Fx_hat[i-1]*x_bar_hat_tmp;
		dcheckX = Fx_hat[i-1]*dx_bar_hat_tmp;
		t_vecX.setConstant(1);
		// cout << "check:" << endl << check << endl << endl;
		for (int j=1; j<= _nSt; j++)
		{
			if(dcheckX[j-1]>0)
				t_vecX[j-1] = checkX[j-1]/dcheckX[j-1];
		}
		t_tmp = t_vecX.minCoeff();
		if (t_tmp < t)
			t = t_tmp;
		
		checkU = fu[i] - (Fu_bar_hat[i-1]*x_bar_hat_tmp + Fu_hat[i]*c_hat_tmp);
		dcheckU = Fu_bar_hat[i-1]*dx_bar_hat_tmp + Fu_hat[i]*dc_hat_tmp;
		t_vecU.setConstant(1);
		for (int j=1; j<=_nInp; j++)
		{
			if(dcheckU[j-1]>0)
				t_vecU[j-1] = checkU[j-1]/dcheckU[j-1];
		}
		t_tmp = t_vecU.minCoeff();
		if (t_tmp < t)
			t = t_tmp;
	}
	// cout << "finished general treatment at the beginning" << endl << endl;
	
	
	// special case for last blocks
	x_bar_hat_tmp = z_hat.template segment<_n+1>((_N-1)*offset2+_m+_n+1);
	dx_bar_hat_tmp = dz_hat.template segment<_n+1>((_N-1)*offset2+_m+_n+1);
	c_hat_tmp = z_hat.template segment<_m+1>(_N*offset2);
	dc_hat_tmp = dz_hat.template segment<_m+1>(_N*offset2);
	checkX = fx[_N-1] - (Fx_hat[_N-1]*x_bar_hat_tmp);
	dcheckX = Fx_hat[_N-1]*dx_bar_hat_tmp;
	t_vecX.setConstant(1);
	for (int j=1; j<= _nSt; j++)
	{
		if (dcheckX[j-1]>0)
			t_vecX[j-1] = checkX[j-1]/dcheckX[j-1];
	}
	t_tmp = t_vecX.minCoeff();
	if (t_tmp < t)
		t = t_tmp;
	
	x_bar_hat_tmp = z_hat.template segment<_n+1>((_pos_omega-1)*offset2+_m+_n+1);
	dx_bar_hat_tmp = dz_hat.template segment<_n+1>((_pos_omega-1)*offset2+_m+_n+1);
	checkTheta = f_xTheta - (F_xTheta_hat*x_bar_hat_tmp + F_theta_hat*c_hat_tmp);
	dcheckTheta = F_xTheta_hat*dx_bar_hat_tmp + F_theta_hat*dc_hat_tmp;
	t_vecTheta.setConstant(1);
	for (int j=1; j<= _nF_xTheta; j++)
	{
		if(dcheckTheta[j-1]>0)
			t_vecTheta[j-1] = checkTheta[j-1]/dcheckTheta[j-1];
	}
	t_tmp = t_vecTheta.minCoeff();
	if (t_tmp < t)
		t = t_tmp;
	
	// return result
	if (t == 1)
		return 1;
	else
		return 0.99*t;	// guarantees strict feasibility	
}



#endif