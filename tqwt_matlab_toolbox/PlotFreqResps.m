function PlotFreqResps(Q, r, J, Fs)
% PlotFreqResps(Q, r, Levels)
% Plot frequency responses of the tunable Q-factor wavelet transform (TQWT)
% with parameters (Q,r).
%
% Use PlotFreqResps(Q, r, Levels, Fs) to specify sampling frequency

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010

check_params(Q,r,J);

figure(gcf)

if nargin < 4
    Fs = 1;
end

beta = 2/(Q+1);
alpha = 1-beta/r;

for j = 1:J
    w2 = (alpha)^(j-1) * pi;
    w1 = (1-beta) * w2;
    w = linspace(w1,w2,70);
    H1j = H1j_fun(w,alpha,beta,j);
    H1j = H1j/max(H1j);
    line(w/pi/2*Fs,H1j);
end
set(gca,'Xtick',(0:0.1:0.5)*Fs)
box off
xlim([0 0.5*Fs])
ylim([0 1.1])
% ylim([0 1.2*max(H1j)])


txt1 = sprintf('FREQUENCY RESPONSES: Q = %1.2f, R = %3.2f', Q,r);
title(txt1)
if nargin < 4
    xlabel('NORMALIZED FREQUENCY (HERZ)')
else
    xlabel('FREQUENCY (HERZ)')
end



% * * * Local Functions * * *  



function H0 = H0_fun(w,alpha,beta)
% H0 = H0_fun(w,alpha,beta)
%
% alpha = 0.8;
% beta = 0.5;
% w = linspace(0,pi,100);
% H0 = H0_fun(w,alpha,beta);
% H1 = H1_fun(w,alpha,beta);
% subplot(2,1,1), plot(w/pi,H0,w/pi,H1)
% subplot(2,1,2), plot(H0.^2 + H1.^2 - 1)

H0 = zeros(size(w));

w = mod(w+pi,2*pi)-pi;

H0(abs(w) <= alpha*pi) = 1;

k = (abs(w) >= (1-beta)*pi) & (abs(w) <= alpha*pi);

H0(k) = theta_fun((w(k)+(beta-1)*pi)/(alpha+beta-1));



function H0 = H0j_fun(w, alpha, beta, j)
% H0 = H0j_fun(w, alpha, beta, j)
% 0 < alpha, beta < 1
% j : non-netative integer

H0 = ones(size(w));
for m = 0:j-1
    H0 = H0 .* H0_fun(w / alpha^m, alpha, beta);
end



function H1 = H1_fun(w, alpha, beta)

% H1 = H1_fun(w, alpha, beta)

H1 = zeros(size(w));

w = mod(w+pi, 2*pi)-pi;

H1(abs(w) >= alpha*pi) = 1;

k = (abs(w) >= (1-beta)*pi) & (abs(w) <= alpha*pi);

H1(k) = theta_fun( (alpha*pi - w(k))/(alpha+beta-1) );



function H1 = H1j_fun(w, alpha, beta, j)

% H1 = H1j_fun(w, alpha, beta, j)
%
% j : non-netative integer

H1 = H0j_fun(w, alpha, beta, j-1) .* H1_fun(w / alpha^(j-1), alpha, beta);



function theta = theta_fun(w)

% theta = theta_fun(w)
% Equation (28) in RADWT paper (Aug 2009, IEEE TSP)

theta = 0.5 * (1 + cos(w)) .* sqrt(2-cos(w));

