function X = uDFT(x)
% X = uDFT(x)
% unitary DFT

N = length(x);
X = fft(x) / sqrt(N);