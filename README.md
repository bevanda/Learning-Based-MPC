lbmpc_ipm
===========

This is an implementation of Learning-Based Model Predictive Control (LBMPC) that uses the sparse LBmpcIPM solver. LBmpcIPM is a free solver written in C++ under the BSD license.

Prerequisites
=============
* LBmpcIPM
* [Eigen3](eigen.tuxfamily.org/)
* CMake

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
