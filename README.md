lbmpc_ipm
===========

This is an implementation of Learning-Based Model Predictive Control (LBMPC) that uses the LBmpcIPM solver. LBmpcIPM is a sparse primal-dual infeasible start interior point method implementation based on Mehrotra's predictor-corrector scheme. This open-source solver is written in C++ and freely available under the BSD licence. The solver is tailored towards QP-LBMPC, and is provided to enable the rapid implementation of LBMPC on other platforms. Another salient feature is that its solving times scales linearly in the prediction horizon.

Aside from the implementation, this repository also contains the source code of LBmpcIPM  and a short documentation on how to use the solver. Further details on the implementation can be found [here](http://control.ee.ethz.ch/index.cgi?page=publications&action=details&id=4168)


Prerequisites
=============

Prior to compiling and running the simulations, make sure the following applications and libraries are available on the computer:

* LBmpcTP.h (template class file of LBmpcIPM)
* [Eigen3](http://eigen.tuxfamily.org/)
* [CMake](http://www.cmake.org/)

Compiling examples
==================

    cd lbmpc_ipm
    mkdir build
    cd build
    cmake ..
    export N_MPC_STEPS=15 # or whatever..
    make

Creating data files
===================

Example from documentation
--------------------------
(in MATLAB):

    >> cd lbmpc_ipm/matlab/example0
    >> Init

Quadrotor example
-----------------
(in MATLAB):

    >> cd lbmpc_ipm/matlab/qr_example
    >> Init

How to run examples
===================

Example from documentation
--------------------------
    cd lbmpc_ipm
    build/bin/example0 matlab/example0/ConstrParam.bin

Quadrotor example
-----------------
	cd lbmpc_ipm
	build/bin/qr_example matlab/qr_example/quad.bin
	
[Back to home](https://bitbucket.org/lbmpc/lbmpc.bitbucket.org/wiki/Home)
