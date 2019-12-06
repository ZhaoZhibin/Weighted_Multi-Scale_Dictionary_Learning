function y = itqwt(w,Q,r,N)
% y = itqwt(w,Q,r,N)
% inverse tunable Q-factor wavelet transform
% See tqwt

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

Y = uDFT(w{J+1});

for j = J:-1:1
    W = uDFT(w{j});
    M = 2*round(alpha^(j-1) * N/2);
    Y = sfb(Y, W, M);
end

y = uDFTinv(Y);

