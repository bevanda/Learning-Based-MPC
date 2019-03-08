%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimal main program for the use with nmpc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function minimal_example
    addpath('./nmpcroutine');
    mpciterations = 10;
    N             = 2;
    T             = 0.1;
    tmeasure      = 0.0;
    xmeasure      = [2.0];
    u0            = zeros(1,N);

    [t, x, u] = nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system, ...
         mpciterations, N, T, tmeasure, xmeasure, u0);
    rmpath('./nmpcroutine');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    cost = x^2+u^2;

function cost = terminalcosts(t, x)
    cost = 0.0;

function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = [];
    ub  = [];

function y = system(t, x, u, T)
    y = x+u;