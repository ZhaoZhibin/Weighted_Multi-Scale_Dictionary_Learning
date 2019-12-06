%% Sparse representation with the TQWT (complex-valued signal)
% Example illustrating sparse signal representation/approximation
% using the tunable Q-factor wavelet transform (TQWT).
% The first part: sparse signal representation (basis pursuit).
% The second part: sparse signal approximation (basis pursuit densoising).
%
% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick,  selesi@poly.edu
% Polytechnic Institute of NYU
% November 2010

%% Miscellaneous

clear
close all

%% Set example parameters

% Uncomment one of the following lines
Example_Num = 1;              % Artificial signal - Low Q-factor
% Example_Num = 2;              % Artificial signal - High Q-factor

switch Example_Num
    case 1
        
        % Set wavelet parameters
        Q = 1;
        r = 3.0;
        L = 9;                 % number of levels
        L1 = 1;
        
        x = test_signal(11);     % Make test signal
        N = length(x);
        t = (0:N-1);            % time axis
        fs = 1;                 % sampling frequency
        xlabel_txt = 'TIME (SAMPLES)';
        A = 2;

        printme = @(str) print('-dpdf',sprintf('figures/sparsity_demo_complex_Ex1_%s',str));
        
    case 2
        
        % Set wavelet parameters
        Q = 3;
        r = 3;
        L = 23;
        L1 = 10;
        
        x = test_signal(12);     % Make test signal
        N = length(x);
        t = (0:N-1);            % time axis
        fs = 1;                 % sampling frequency
        xlabel_txt = 'TIME (SAMPLES)';
        A = 2;

        printme = @(str) print('-dpdf',sprintf('figures/sparsity_demo_complex_Ex2_%s',str));

end

beta = 2/(Q+1);
alpha = 1-beta/r;


%% Display signal

figure(1), clf
subplot(2,1,1)
plot(t, real(x), t, imag(x),'r')
legend('REAL PART','IMAGINARY PART')
title('TEST SIGNAL')
box off
ylim([-A A])
xlim([t(1) t(end)])
xlabel(xlabel_txt)

% Verify perfect reconstruction

w = tqwt_radix2(x, Q, r, L);
y = itqwt_radix2(w, Q, r, N);

fprintf('Reconstruction error = %4.3e\n', max(abs(y-x)))



%% Plot wavelet at final level and low-pass scaling function

wlets = ComputeWavelets(N,Q,r,L,'radix2');      % Compute wavelets

figure(3), clf
subplot(2,1,1)
plot(t, wlets{L})
title(sprintf('WAVELET AT LEVEL %d',L))
xlim([0 N/fs])
ylim(1.2*(max(abs(wlets{L})))*[-1 1])
box off
subplot(2,1,2)
plot(t, wlets{L+1})
title(sprintf('LOW-PASS SCALING FUNCTION AT LEVEL %d',L))
xlim([0 N/fs])
box off
xlabel(xlabel_txt)


%% Plot subbands

figure(4), clf
PlotSubbands(x,w,Q,r,L1,L,fs,'E');
title('SUBBANDS')
orient tall
printme('subbands')


%% Compute energy in each subband
% It can be useful to know how the energy of a signal is distributed
% across the subbands. We compute the energy in each subband and display
% using a bar graph. Because the transform has the Parseval property the
% distribution of the energy across the subbands reflects the frequency
% content of the signal.

figure(5)
clf
subplot(3,1,1:2)
e = PlotEnergy(w);
printme('energy')


%% Sparse wavelet representation (Basis Pursuit)
% Sparse signal representation with perfect reconstruction (Basis Pursuit).
% Use a variant of SALSA to minimize l1-norm of wavelet coefficients
% providing exact reconstruction.

now = ComputeNow(N,Q,r,L,'radix2');

lambda = now;       % Regularization parameters
mu = 1.0;           % SALSA parameter
Nit = 200;          % Number of iterations

[w2, costfn] = tqwt_bp(x, Q, r, L, lambda, mu, Nit);
y = itqwt_radix2(w2, Q, r, N);

err = x - y;

rel_err = sqrt(mean(err.^2))/sqrt(mean(x.^2));

fprintf('Basis pursuit (BP) relative RMS reconstruction error = %4.3e\n',rel_err);

%% Compute cost function

cost = 0;
for j = 1:L+1
    cost = cost + lambda(j)*sum(abs(w2{j}));
end

fprintf('BP objective function: %e\n', cost);

%% Display cost function versus iteration

figure(6), clf
it1 = 10;
plot(it1:Nit, costfn(it1:Nit));
xlim([0 Nit])
title('SALSA COST FUNCTION (BASIS PURSUIT)')
xlabel('ITERATION')
box off
printme('cost_function_bp')

% check consistency between 'cost' and final value of 'costfn'
fprintf('BP costfn(end) = %d\n', costfn(end))


%% Plot sparse subbands

figure(7), clf
PlotSubbands(y,w2,Q,r,L1,L,fs,'E');
title('SPARSE SUBBANDS (BASIS PURSUIT)')
orient tall
printme('subbands_sparse_bp')

%% Reconstruction error

figure(8), clf
subplot(2,1,1)
plot(t, real(x), t, imag(x))
xlim([t(1) t(end)])
ylim([-A A])
box off
title('TEST SIGNAL')

subplot(2,1,2)
plot(t, abs(x-y))
xlim([t(1) t(end)])
box off
title('RECONSTRUCTION ERROR (BASIS PURSUIT)')
xlabel(xlabel_txt)

printme('recon_error_bp')

%% Compute energy in each subband

figure(9)
clf
subplot(3,1,1:2)
e2 = PlotEnergy(w2);
title('DISTRIBUTION OF SIGNAL ENERGY (BASIS PURSUIT)')
printme('energy_bp')



%% Sparse wavelet approximation (Basis Pursuit Denoising)
% Sparse signal representation with l1-norm regularization (Basis Pursuit Denoising).
% Use SALSA to minimize function: sum((x-invTQWT(w)).^2) + sum(abs((lambda.*w))).

lambda = 0.002*now;         % Regularizaton parameter
mu = 0.0004;                % SALSA parameter
Nit = 200;                  % Number of iterations

[w3, costfn] = tqwt_bpd(x, Q, r, L, lambda, mu, Nit);
y = itqwt_radix2(w3, Q, r, N);

err = x - y;

rel_err = sqrt(mean(err.^2))/sqrt(mean(x.^2));

fprintf('Basis pursuit denoising (BPD) relative RMS reconstruction error = %4.3e\n',rel_err);


%% Compute cost function

cost = sum(abs(x - y).^2);
for j = 1:L+1
    cost = cost + lambda(j)*sum(abs(w3{j}));
end

fprintf('BPD objective function: %e\n', cost);

%% Display cost function versus

figure(11), clf
it1 = 10;
plot(it1:Nit, costfn(it1:Nit));
xlim([0 Nit])
title('SALSA COST FUNCTION (BASIS PURSUIT DENOISING)')
xlabel('ITERATION')
box off
printme('cost_function_bpd')

% check consistency between 'cost' and final value of 'costfn'
fprintf('BPD costfn(end) = %d\n', costfn(end))


%% Plot sparse subbands

figure(10), clf
PlotSubbands(y,w3,Q,r,L1,L,fs,'E');
title('SPARSE SUBBANDS (BASIS PURSUIT DENOISING)')
orient tall
printme('subbands_sparse_bpd')

%% Reconstruction error

figure(11), clf
subplot(2,1,1)
plot(t, real(x), t, imag(x))
xlim([t(1) t(end)])
ylim([-A A])
box off
title('TEST SIGNAL')

subplot(2,1,2)
plot(t, abs(x-y))
xlim([t(1) t(end)])
box off
% ylim([-A A]*0.1)
title('RECONSTRUCTION ERROR (BASIS PURSUIT DENOISING)')
xlabel(xlabel_txt)

printme('recon_error_bpd')

%% Write information to file

file_name = sprintf('figures/sparsity_demo_complex_Ex%d_info.txt',Example_Num);
fid = fopen(file_name,'w');
fprintf(fid,'TRANSFORM PARAMETERS:\n');
fprintf(fid,'\t Transform: tqwt_radix2.m\n');
fprintf(fid,'\t Q = %4.2f\n\t r = %4.2f\n\t levels = %d\n\n',Q,r,L);
fprintf(fid,'SALSA PARAMETERS: \n');
fprintf(fid,'\t mu = %.2e\n\t iterations = %d\n', mu, Nit);
fclose(fid);
