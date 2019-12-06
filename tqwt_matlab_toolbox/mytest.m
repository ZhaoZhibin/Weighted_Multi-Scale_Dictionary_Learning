Q = 9; r = 3; J = 10; % TQWT parameters
figure(1), clf
N = 256;
PlotWavelets(N, Q, r, 1, J, 'radix2');


x = test_signal(4); % Make test signal
N = length(x); % Length of test signal
w = tqwt_radix2(x,Q,r,J); % TQWT
y = itqwt_radix2(w,Q,r,N); % Inverse TQWT

% w = tqwt(x,Q,r,J); % TQWT
% y = itqwt(w,Q,r,N); % Inverse TQWT

recon_err = max(abs(x - y)) % Reconstruction error