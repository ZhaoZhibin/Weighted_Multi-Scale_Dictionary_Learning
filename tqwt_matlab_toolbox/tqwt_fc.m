function fc = tqwt_fc(Q, r, J, fs)
% fc = tqwt_fc(Q, r, J, fs)
% fc - center frequencies of each subband
% fc(k) is center frequency of subband k (k = 1..J)

check_params(Q,r,J);

beta = 2/(Q+1);
alpha = 1-beta/r;
fc = [0.5 alpha.^(2:J) * (2-beta)/(4*alpha)] * fs;
