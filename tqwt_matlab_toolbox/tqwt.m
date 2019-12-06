function w = tqwt(x,Q,r,J)
% w = tqwt(x,Q,r,J)
% Tunable Q-factor wavelet transform (TQWT)
% INPUT
%   x - input signal (even length)
%   Q - Q-factor
%   r - oversampling rate (redundancy)
%   J - number of levels
% OUTPUT
%   w - wavelet coefficients
% NOTES
%   w{j} is subband j for j = 1,...,J+1
%   w{1} is the highest-frequency subband signal
%   w{J+1} is the low-pass subband signal
%
% % Example (verify perfect reconstruction)
% Q = 4; r = 3; J = 3;     % parameters
% N = 200;                 % signal length
% x = rand(1,N);           % test signal
% w = tqwt(x,Q,r,J);       % wavelet transform
% y = itqwt(w,Q,r,N);      % inverse wavelet transform
% max(abs(x - y))          % reconstruction error

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010

check_params(Q,r,J);

beta = 2/(Q+1);
alpha = 1-beta/r;
N = length(x);

Jmax = floor(log(beta*N/8)/log(1/alpha));
if J > Jmax
    J = Jmax;
    fprintf('Note: too many levels specified.\n')
    if Jmax > 0
        fprintf('Reduce levels to %d\n',Jmax);
    else
        fprintf('Increase signal length.\n');
    end
    w = [];
    return
end

X = uDFT(x);
w = cell(1, J+1);

% J stages:
for j = 1:J
    N0 = 2*round(alpha^j * N/2);
    N1 = 2*round(beta * alpha^(j-1) * N/2);
    [X, W] = afb(X, N0, N1);
    w{j} = uDFTinv(W);
end

w{J+1} = uDFTinv(X);

