function Jmax = MaxLevels(Q,r,N)
% Jmax = MaxLevels(Q,r,N)
% Maximum number of levels of TQWT with parameters (Q,r) for an N-point signal.

beta = 2/(Q+1);
alpha = 1-beta/r;

Jmax = floor(log(beta*N/8)/log(1/alpha));