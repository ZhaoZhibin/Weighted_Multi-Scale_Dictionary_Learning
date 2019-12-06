function now = ComputeNow_radix2(N,Q,r,J)

% Compute norms of wavelets (now) for radix-2 TQWT
%
% now = ComputeNow_radix2(N,Q,r,J)
%
% now(j) is the L2 norm of the wavelet at level j, for j = 1,...,J+1
% for an N-point radix-2 TQWT with parameters (Q,r).
%
% Note: N must be a power of 2.

% This is an mfile version of the the mex function.
% The mex function will be faster if it is available and
% functional for your Matlab version and operating system.

% If no mex function is available, then 'now' can be computed
% using the mfile ComputeNow.m:

now = ComputeNow(N,Q,r,J,'radix2');


