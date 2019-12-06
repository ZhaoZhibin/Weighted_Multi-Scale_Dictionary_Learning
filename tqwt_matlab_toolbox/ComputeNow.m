function now = ComputeNow(N,Q,r,J,radix_flag)

% Compute norms of wavelets (now)
%
% now = ComputeNow(N,Q,r,J)
% now = ComputeNow(N,Q,r,J,'radix_2')
%
% now(j) is the L2 norm of the wavelet at level j, for j = 1,...,J+1
% for an N-point TQWT with parameters (Q,r).
%
% The norms of the wavelets are needed when performing wavelet-domain
% thresholding.
%
% Note: for 'radix2' option, N must be a power of 2.


% This program produces the same 'now' vector as 'ComputeWavelets',
% however, this program is more efficient because it does not use
% the FFT. Because of Parseval's theorem, the norm of the wavelet
% can be computed in the frequency domain.

% Check that the two functions produce the same 'now' vector:
if 0
    Q = 3; r = 3; J = 10; N = 300;
    [wlets, now1] = ComputeWavelets(N, Q, r, J);
    [now2] = ComputeNow(N, Q, r, J);
    now1 - now2
    
    Q = 3; r = 3; J = 10; N = 256;
    [wlets, now1] = ComputeWavelets(N, Q, r, J, 'radix2');
    [now2] = ComputeNow(N, Q, r, J, 'radix2');
    now1 - now2
end

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010



if nargin == 5
    if strcmp(radix_flag, 'radix2')
        R2 = true;
    else
        disp('invalid string')
        now = [];
        return
    end
else
    R2 = false;
end

if R2 == true
    if log2(N) ~= round(log2(N))
        disp('N must be a power of 2 for radix-2 option for computing norm of wavelets')
        disp('(otherwise, not all wavelets in each subband has equal norm).')
        now = [];
        return
    end
end



beta = 2/(Q+1);
alpha = 1-beta/r;


% Create all-zero wavelet coefficients
w = cell(1, J+1);
for j = 1:J
    N0 = 2*round(alpha^j * N/2);
    N1 = 2*round(beta * alpha^(j-1) * N/2);
    if R2 == true
        w{j} = zeros(1,next(N1));
    else
        w{j} = zeros(1,N1);
    end
end
if R2 == true
    w{J+1} = zeros(1,next(N0));
else
    w{J+1} = zeros(1,N0);
end

wz = w;

now = zeros(1,J+1);


for i = 1:J+1    
    
    w = wz;
    
    %     w{i}(1) = 1;
    M = length(w{i});
    w{i}(1:M) = 1/sqrt(M);         % define directly in uDFT, then no uDFT is needed.
    
    %     Y = uDFT(w{J+1});
    Y = w{J+1};
    
    if R2 == true
        M = 2*round(alpha^J * N/2);     % *
        Y = lps(Y, M);                  % *
    end
    
    for j = J:-1:1
        %         W = uDFT(w{j});
        W = w{j};
        if R2 == true
            N1 = 2*round(beta * alpha^(j-1) * N/2);     % *
            W = lps(W, N1);                             % *
        end
        
        M = 2*round(alpha^(j-1) * N/2);
        Y = sfb(Y, W, M);
    end
    
    now(i) = sqrt(sum(abs(Y).^2));       % l2 norm of wavelet
    
end

