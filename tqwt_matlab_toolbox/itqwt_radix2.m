function y = itqwt_radix2(w,Q,r,L)
% y = itqwt_radix2(w,Q,r,N)
% inverse radix-2 TQWT
% See tqwt_radix2

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010

check_params(Q,r);

beta = 2/(Q+1);
alpha = 1-beta/r;
J = length(w)-1;

N = next(L);           % *

Y = uDFT(w{J+1});

M = 2*round(alpha^J * N/2);     % *
Y = lps(Y, M);                  % *

for j = J:-1:1
    W = uDFT(w{j});
    N1 = 2*round(beta * alpha^(j-1) * N/2);     % *
    W = lps(W, N1);                             % *
    M = 2*round(alpha^(j-1) * N/2);
    Y = sfb(Y, W, M);
end

y = uDFTinv(Y);
y = y(1:L);                     % *

% * : denotes modification for radix-2 case

