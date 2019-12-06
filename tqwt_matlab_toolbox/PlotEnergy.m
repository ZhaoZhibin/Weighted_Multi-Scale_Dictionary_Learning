function e = PlotEnergy(w, Q, r, fs)
% e = PlotEnergy(w)
% Plot distribution of signal energy over subbands
%
% Use PlotEnergy(w, Q, r, fs) to plot energy versus frequency 
% (the low-pass band is excluded from the plot)
%
% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
%
% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010


J = length(w)-1;

e = zeros(1,J+1);
for j = 1:J+1
    e(j) = sum(abs(w{j}).^2);
end
figure(gcf)


if nargin == 1
    
    b1 = bar(100*e/sum(e));
    set(b1,'facecolor','w')
    xlim([0 J+2])
    xlabel('SUBBAND')
    ylabel('SUBBAND ENERGY (% OF TOTAL)')
    
    
else

    beta = 2/(Q+1);
    alpha = 1-beta/r;
    fc = [alpha.^(J:-1:2) * (2-beta)/(4*alpha) 0.5] * fs;       % center frequencies

    fb = [alpha.^((J:-1:1)+0.5) * (2-beta)/(4*alpha) 0.5] * fs;       % bandedges frequencies

    % make frequency-axis tick marks
    m = ceil(log2(fc(end)/fc(1)));
    f_tick = fs * 2.^(-m:1:-1);

    ee = [e(J:-1:1) e(1)];
    [xx,yy] = stairs(fb, 100*ee/sum(e));
    
    % semilogx(xx, yy);
    fill(([xx' xx(end) xx(1)]), [yy' 0 0], [1 1 1]*0.7)
    set(gca,'Xscale','log')
    
    xlim([fb(1) fb(end)])
    set(gca,'Xtick', f_tick)
    xlabel('FREQUENCY (HZ)')
    grid
    ylabel('ENERGY (% OF TOTAL)')
    
end

title('DISTRIBUTION OF SIGNAL ENERGY')
box off

