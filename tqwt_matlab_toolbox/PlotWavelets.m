function wlets = PlotWavelets(N,Q,r,J1,J2,radix_flag)
% wlets = PlotWavelets(N,Q,r,J1,J2)
% Plots TQWT wavelets from scale J1 to J2 for an N-point signal.
% Use PlotWavelets(N,Q,r,J1,J2,'radix2') for radix-2 TQWT wavelets.
%
% % Example:
% N = 2^9, Q = 3.5, r = 3, J = 20;
% PlotWavelets(N,Q,r,5,J,'radix2');
% % PlotWavelets(N,Q,r,5,J);

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010

check_params(Q,r,J1);

beta = 2/(Q+1);
alpha = 1-beta/r;

if nargin  == 6
    wlets = ComputeWavelets(N,Q,r,J2,radix_flag);
else
    wlets = ComputeWavelets(N,Q,r,J2);
end

if isempty(wlets)
    return
end

figure(gcf)
clf
for j = J1:J2
    y = wlets{j};
    M = 2.5*max(abs(y));
    plot((1:N),y/M-j)            % Plot wavelet
    hold on
    
    if 0
        % To show (approx) support of wavelet
        BW = alpha^(j-1) * beta / 4;    % bandwidth (cycles/sample)
        dur = 2 / BW;                   % approximate duration
        [mm,km] = max(abs(y));
        plot(km+[-1 1]*dur/2,[1 1]*mm/M*0.2-j,'r')
    end
end
ylim([-J2-1 -J1+1])

hold off
xlim([0 N])
xlabel('TIME (SAMPLES)')

set(gca,'Ytick',-J2:-J1,'YtickLabel', J2:-1:J1)
ylabel('SUBBAND')
title_string = sprintf('WAVELETS: SUBBANDS %d THROUGH %d',J1,J2);
title(title_string)

txt = sprintf('Q = %.2f, r = %.2f',Q,r);
% text(0,-0.05,txt,'units','normalized');

tb = annotation('textbox',[.0 .0 .6 .05]);
set(tb,...
    'string',txt,...
    'verticalalignment','bottom',...
    'edgecolor','none')

%     'horizontalalignment','center',...



