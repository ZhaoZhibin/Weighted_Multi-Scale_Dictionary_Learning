function [x1,x2,w1,w2,costfn] = dualQ(x,Q1,r1,J1,Q2,r2,J2,lam1,lam2,mu,Nit,plot_flag)
% [x1,x2,w1,w2] = dualQ(x,Q1,r1,J1,Q2,r2,J2,lam1,lam2,mu,Nit)
% Resonance signal decomposition using two Q-factors.
% This function minimizes the cost function:
%   lam1||w1||_1 + lam2||w2||_2  subject to x = TQWT(w1) + TQWT(w2)
% (see guide for clarification).
%
% INPUT
%   x - input signal signal
%   Q1, r1, J1, Q1, r1, J2 - TQWT parameters
%   lam1, lam2 - regularization parameters
%   mu - SALSA parameter
%   Nit - Number of iterations
% OUTPUT
%   x1, x2 - components
%   w1, w2 - transform coefficients of components
%
% Use [x1,x2,w1,w2,costfn] = dualQ(...) to return cost function.
%
% Use [...] = dualQ(...,'plots') to plot progress of algorithm (this
% option is not available in the mex version).
%
% This function uses the radix-2 TQWT.  If length(x) is not a
% power of 2, then x will be zero-padded to next power of 2.
%
% Algorithm: SALSA-type algorithm for basis pursuit with TQWT.

% References: 
%
% Resonance-Based Signal Decomposition: A New Sparsity-Enabled Signal Analysis Method. 
% Signal Processing, 2010, doi:10.1016/j.sigpro.2010.10.018
% I. W. Selesnick. 
%
% To minimize the cost function, this program uses SALSA. SALSA is given in:
% M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo.
% Fast image recovery using variable splitting and constrained optimization.
% IEEE Trans. Image Process., 19(9):2345-2356, September 2010.

check_params(Q1,r1,J1);
check_params(Q2,r2,J2);

% By default do not compute cost function (to reduce computation)
if nargout > 4
    COST = true;
    costfn = zeros(1,Nit);     % cost function
else
    COST = false;
    costfn = [];
end

GOPLOTS = false;
if nargin == 12
    if strcmp(plot_flag,'plots')
        GOPLOTS = true;
    end    
end

L = length(x);
N = next(L);
if L < N
    x = zeropad(x,N);
end

x = x(:).';         % row vector

% Initialize:
w1 = tqwt_radix2(x,Q1,r1,J1);
w2 = tqwt_radix2(x,Q2,r2,J2);
d1 = tqwt_radix2(zeros(size(x)),Q1,r1,J1);
d2 = tqwt_radix2(zeros(size(x)),Q2,r2,J2);

T1 = lam1/(2*mu);
T2 = lam2/(2*mu);

u1 = cell(1,J1+1);
u2 = cell(1,J2+1);

N = length(x);
A = 1.1*max(abs(x));

for k = 1:Nit
    
    for j = 1:J1+1
        u1{j} = soft(w1{j} + d1{j}, T1(j)) - d1{j};
    end
    for j = 1:J2+1
        u2{j} = soft(w2{j} + d2{j}, T2(j)) - d2{j};
    end
    
    c = x - itqwt_radix2(u1,Q1,r1,N) - itqwt_radix2(u2,Q2,r2,N);
    c = 0.5 * c;
    
    d1 = tqwt_radix2(c,Q1,r1,J1);
    d2 = tqwt_radix2(c,Q2,r2,J2);
    
    for j = 1:J1+1
        w1{j} = d1{j} + u1{j};
    end
    
    for j = 1:J2+1
        w2{j} = d2{j} + u2{j};
    end
    
    if COST
        % Note: due to equality constraint, residual is zero.
        % So the residual is not included in the cost function.
        costfn(k) = 0;
        for j = 1:J1+1
            costfn(k) = costfn(k) + lam1(j)*sum(abs(w1{j}));
        end
        for j = 1:J2+1
            costfn(k) = costfn(k) + lam2(j)*sum(abs(w2{j}));
        end
    end
    
    if GOPLOTS
        x1 = itqwt_radix2(w1,Q1,r1,N);
        x2 = itqwt_radix2(w2,Q2,r2,N);
        res = x - x1 - x2;

        figure(gcf)
        clf
        subplot(3,1,1)
        plot(x1)
        xlim([0 N])
        ylim([-A A])
        title({sprintf('ITERATION %d',k),'COMPONENT 1'})
        box off
        subplot(3,1,2)
        plot(x2)
        xlim([0 N])
        ylim([-A A])
        box off
        title('COMPONENT 2')
        subplot(3,1,3)
        plot(res)
        xlim([0 N])
        ylim([-A A])
        title('RESIDUAL')
        box off
        drawnow
    end
    
end

x1 = itqwt_radix2(w1,Q1,r1,L);
x2 = itqwt_radix2(w2,Q2,r2,L);


% --------------- local function: soft ---------------


function y = soft(x,T)
% Soft-threshold function
% y = soft_fun(x,T)
% x : input data
% T : threshold

if isreal(x)
    y = zeros(size(x));
    k = (x < -T);
    y(k) = x(k) + T;
    k = (x > T);
    y(k) = x(k) - T;
else
    % following alternative definition works for real and complex data:
    y = max(abs(x)-T,0);
    y = y./(y+T) .* x;
end



