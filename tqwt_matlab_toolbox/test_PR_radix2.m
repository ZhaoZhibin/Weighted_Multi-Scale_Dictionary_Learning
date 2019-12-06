%% Radix-2 Tunable Q-factor Wavelet Transform: Verify PR
% Verify perfect reconstruction (PR) property of radix-2 TQWT

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick,  selesi@poly.edu
% Polytechnic Institute of NYU
% November 2010

clear

for k = 1:2

    if k == 1
        Q = 4; r = 3; J = 10;   % High Q-factor wavelet transform
    elseif k == 2
        Q = 1; r = 3; J = 5;    % Low Q-factor wavelet transform
    end

    fprintf('\n')
    fprintf('Q = %3.2f, r = %3.2f\n',Q,r)
    
    beta = 2/(Q+1);
    alpha = 1-beta/r;
    I = sqrt(-1);

    % Verify PR
    for N = 2.^(7:10)                   % Verify PR for various lengths
        x = rand(1,N) + I*rand(1,N);    % Make test signal (complex-valued)
        J = floor(log2(beta*N/8)/log2(1/alpha));    % number of levels
        w = tqwt_radix2(x,Q,r,J);       % TQWT
        y = itqwt_radix2(w,Q,r,N);      % Inverse TQWT
        recon_err = max(abs(x - y));    % Reconstruction error

        fprintf('N = %4d, J = %3d: tqwt/itqwt (radix2) recon error = %e\n',N,J,recon_err)
    end

    for N = 600:2:620                   % Verify PR for various lengths
        x = rand(1,N) + I*rand(1,N);    % Make test signal (complex-valued)
        J = floor(log2(beta*N/8)/log2(1/alpha));    % number of levels
        w = tqwt_radix2(x,Q,r,J);       % TQWT
        y = itqwt_radix2(w,Q,r,N);      % Inverse TQWT
        recon_err = max(abs(x - y));    % Reconstruction error

        fprintf('N = %4d, J = %3d: tqwt/itqwt (radix2) recon error = %e\n',N,J,recon_err)
    end


end


