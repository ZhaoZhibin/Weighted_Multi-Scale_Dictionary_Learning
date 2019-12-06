function x = uDFTinv(X)
% x = uDFTinv(X)
% inverse unitary DFT

N = length(X);
x = sqrt(N) * ifft(X);


