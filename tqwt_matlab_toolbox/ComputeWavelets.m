function [wlets, now] = ComputeWavelets(N,Q,r,J,radix_flag)
% wlets = ComputeWavelets(N,Q,r,J)
% For the tunable Q-factor wavelet transform (TQWT) this function
% computes the N-point wavelets for a J-level transform.
% 
% Use ComputeWavelets(N,Q,r,J,'radix2') for radix-2 TQWT
%
% [wlets, now] = ComputeWavelets(...) also returns a vector of
% the norms of the wavelets for each subband.
%
% % Example:
% Q = 1; r = 3; J = 10; N = 2^9;   % or
% Q = 4; r = 3; J = 18; N = 2^9;
% wlets = ComputeWavelets(N,Q,r,J,'radix2');
% figure(1), clf,
% for j = 1:J+1, line(1:N,wlets{j}/max(abs(wlets{j}))/2+j); end; xlim([0 N])

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010

if nargin == 5
    if strcmp(radix_flag, 'radix2')
        xform = @tqwt_radix2;
        inv_xform = @itqwt_radix2;
        C = N/next(N);              % reduce according to amount of zero-padding
    else
        disp('invalid string')
        wlets = [];
        return
    end
else
    xform = @tqwt;
    inv_xform = @itqwt;
    C = 1;
end

z = zeros(1,N);                     % All-zero signal

wz = xform(z,Q,r,J);                % All-zero wavelet coefficients

if isempty(wz)
    wlets = [];
    return
end

wlets = cell(1,J+1);

for j = 1:J+1
    w = wz;                         % Set w to all-zero coefficients
    m = round(C*length(w{j})/2)+1;  % m: index of coefficient corresponding to midpoint of signal   
    w{j}(m) = 1;                    % Set single wavelet coeff to 1
    y = inv_xform(w,Q,r,N);         % Inverse TQWT
    wlets{j} = y;
end

now = zeros(1,J+1);
for j = 1:J+1
    now(j) = sqrt(sum(abs(wlets{j}).^2));       % l2 norm of wavelet
end


