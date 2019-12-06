function [w, costfn, y] = tqwt_bpd(x, Q, r, J, lam, mu, Nit)
% w = tqwt_bpd(x, Q, r, J, lam, mu, Nit)
% Sparse signal approximation (basis pursuit denoising) using the TQWT.
% Find w to minimize ||x-invTQWT(w)||_2^2 + ||lam.*w||_1 so as to find a
% sparse approximation of signal x.
%
% See sparsity_demo.m
%
% INPUT
%   x - signal
%   Wavelet parameters:
%     Q - Q-factor
%     r - oversampling rate (redundancy)
%     J - number of levels
%   SALSA parameters:
%     lam - regularization parameter (vector)
%     mu  - SALSA parameter
%     Nit - Number of iterations
%
% OUTPUT
%   w - wavelet coefficients
%
% Use [w, costfn] = tqwt_bpd(...) to obtain the cost function
% for each iteration (to see convergence behaviour).  Computing the 
% cost function increases run-time.
%
% Notes: This function uses the radix-2 TQWT.  If length(x) is not a
% power of 2, then x will be zero-padded to next power of 2.

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick,  selesi@poly.edu
% Polytechnic Institute of NYU
% November 2010

% To minimize the cost function, this program uses SALSA.
% SALSA is described in the paper:
% M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo.
% Fast image recovery using variable splitting and constrained optimization.
% IEEE Trans. Image Process., 19(9):2345-2356, September 2010.

check_params(Q,r,J);

% By default do not compute cost function (to reduce computation)
if nargout > 1
    COST = true;
    costfn = zeros(1,Nit);     % cost function
else
    COST = false;
    costfn = [];
end

L = length(x);
N = next(L);
if L < N
    x = zeropad(x,N);
end

% Initialize:
w = tqwt_radix2(x,Q,r,J);
d = tqwt_radix2(zeros(size(x)),Q,r,J);

T = lam/(2*mu);

u = cell(1,J+1);

C = 1/(mu+1);

% SALSA iterations:

for k = 1:Nit
    
    for j = 1:J+1
        u{j} = soft(w{j}+d{j}, T(j)) - d{j};
    end
    
    d = tqwt_radix2(C*x - C*itqwt_radix2(u,Q,r,N),Q,r,J);
        
    for j = 1:J+1
        w{j} = d{j} + u{j};
    end
    
    % If cost function is to be computed...
    if COST
        res = x - itqwt_radix2(w,Q,r,N);              % residual
        costfn(k) = sum(abs(res).^2);
        for j = 1:J+1
            costfn(k) = costfn(k) + lam(j)*sum(abs(w{j}));
        end
    end
end

if nargout > 2
    y = itqwt_radix2(w,Q,r,N);          % for debugging...
end

% --- local functions ----

function y = soft(x,T)
% Soft-threshold function
% y = soft_fun(x,T)
% x : input data
% T : threshold

if isreal(x)
    y = zeros(size(x));
    k = (x < -T);
    y(k) = x(k) + T;
    k = (x > T);
    y(k) = x(k) - T;
else
    % following alternative definition works for real and complex data:
    y = max(abs(x)-T,0);
    y = y./(y+T) .* x;
end


