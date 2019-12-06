function [V0, V1] = afb(X, N0, N1)
% [V0, V1] = afb(X, N0, N1)
% afb: analysis filter bank
% Converts a vector X into two vectors:
% V0 of length N0, V1 of length N1
%
% Need:
%   N0, N1, both even integers
%   2 <= N0 <= length(x)
%   2 <= N1 <= length(x)
%   N0 + N1 > length(x)
%
% % Example (verify perfect reconstruction)
% X = rand(1,20);
% [V0, V1] = afb(X, 16, 12);
% Y = sfb(V0, V1, 20);
% max(abs(X - Y))

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU

% November 2010

X = X(:).';                                     % X is row vector
N = length(X);

P = (N-N1)/2;
T = (N0+N1-N)/2 - 1;
S = (N-N0)/2;

% transition-band function
v = (1:T)/(T+1)*pi;
trans = (1+cos(v)) .* sqrt(2-cos(v))/2;

% Add 1 to indices because Matlab indexing starts at 1 (not 0)

% low-pass subband
V0 = zeros(1,N0);
V0(0+1) = X(0+1);                               % dc term
V0((1:P)+1) = X((1:P)+1);                       % pass-band (pos freq)
V0(P+(1:T)+1) = X(P+(1:T)+1).*trans;            % trans-band (pos freq)
V0(N0/2+1) = 0;                                 % Nyquist freq
V0(N0-P-(1:T)+1) = X(N-P-(1:T)+1).*trans;       % trans-band (neg freq)
V0(N0-(1:P)+1) = X(N-(1:P)+1);                  % pass-band (neg freq)

% high-pass subband
V1 = zeros(1,N1);
V1(0+1) = 0;                                    % dc term
V1((1:T)+1) = X(P+(1:T)+1).*trans(T:-1:1);      % trans-band (pos freq)
V1(T+(1:S)+1) = X(P+T+(1:S)+1);                 % pass-band (pos freq)
if rem(N,2) == 0
    V1(N1/2+1) = X(N/2+1);                      % Nyquist freq (if N even)
end
V1(N1-T-(1:S)+1) = X(N-P-T-(1:S)+1);            % pass-band (neg freq)
V1(N1-(1:T)+1) = X(N-P-(1:T)+1).*trans(T:-1:1); % trans-band (neg freq)




