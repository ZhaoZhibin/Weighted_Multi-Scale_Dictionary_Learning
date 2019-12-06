function Y = sfb(V0, V1, N)
% Y = sfb(V0, V1, N)
% sfb: synthesis filter bank

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU

% November, 2010


N0 = length(V0);
N1 = length(V1);

S = (N-N0)/2;
P = (N-N1)/2;
T = (N0+N1-N)/2 - 1;

% transition-band function
v = (1:T)/(T+1)*pi;
trans = (1+cos(v)) .* sqrt(2-cos(v))/2;

% Add 1 to indices because Matlab indexing starts at 1 (not 0)

% low-pass subband
Y0 = zeros(1,N);
Y0(0+1) = V0(0+1);                                  % dc term
Y0((1:P)+1) = V0((1:P)+1);                          % pass-band (pos freq)
Y0(P+(1:T)+1) = V0(P+(1:T)+1).*trans;               % trans-band (pos freq)
Y0(P+T+(1:S)+1) = 0;                                % stop-band (pos freq)
if rem(N,2) == 0
    Y0(N/2+1) = 0;                                  % Nyquist freq (if N even)
end
Y0(N-P-T-(1:S)+1) = 0;                              % stop-band (neg freq)
Y0(N-P-(1:T)+1) = V0(N0-P-(1:T)+1).*trans;          % trans-band (neg freq)
Y0(N-(1:P)+1) = V0(N0-(1:P)+1);                     % pass-band (neg freq)

% high-pass subband
Y1 = zeros(1,N);
Y1(0+1) = 0;                                        % dc term
Y1((1:P)+1) = 0;                                    % stop-band (pos freq)
Y1(P+(1:T)+1) = V1((1:T)+1).*trans(T:-1:1);         % trans-band (pos freq)
Y1(P+T+(1:S)+1) = V1(T+(1:S)+1);                    % pass-band (pos freq)
if rem(N,2) == 0
    Y1(N/2+1) = V1(N1/2+1);                         % Nyquist freq (if N even)
end
Y1(N-P-T-(1:S)+1) = V1(N1-T-(1:S)+1);               % pass-band (neg freq)
Y1(N-P-(1:T)+1) = V1(N1-(1:T)+1).*trans(T:-1:1);    % trans-band (neg freq)
Y1(N-(1:P)+1) = 0;                                  % stop-band (neg freq)

Y = Y0 + Y1;

