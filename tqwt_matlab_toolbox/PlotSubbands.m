function PlotSubbands(x,w,Q,r,J1,J2,fs,eflag,stem_flag,param_flag)
% PlotSubbands(x,w,Q,r,J1,J2,fs)
% Plot subbands from level J1 to J2 of the TQWT.
% w should be the output of tqwt(x,Q,r,J).
% Need:
%   1 <= J1 < J2 <= J+1 where J is total number of levels.
%
% PlotSubbands(x,w,Q,r,J1,J2,fs,'energy') - displays the
% subband energy as a percentage of the total energy.
%
% PlotSubbands(x,w,Q,r,J1,J2,fs,[],'stem') - uses stem plot style

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU
% November, 2010

N = length(x);

if isreal(x)
    plot((0:N-1)/fs,x/(2*max(abs(x))),'b')
else
    plot((0:N-1)/fs, real(x)/(2*max(abs(x))), 'black')
    line((0:N-1)/fs, imag(x)/(2*max(abs(x))), 'color','blue')
end
% hold on
J = length(w)-1;

ENERGY = false;
if nargin > 7
    if ~isempty(eflag)
        ENERGY = true;
        e = zeros(1,J+1);
        for j = 1:J+1
            e(j) = sum(abs(w{j}).^2);
        end
    end
end

STEM = false;
if nargin > 8
    if strcmp(stem_flag,'stem')
        STEM = true;
    end
end

PARAM = true;
if nargin > 9
    PARAM = false;
end

Mj = zeros(1,J+1);
for j = J1:min([J2 J])
    Mj(j) = max(abs(w{j}));
end
Mj = 0.5 / max(Mj);


% Determine if w is from radix-2 version of TQWT.
% (Yes, if all subbands are of length power 2)
radix2_flag = true;
for j = 1:J+1
    Nj = length(w{j});
    if abs(round(log2(Nj)) - log2(Nj)) > eps
        radix2_flag = false;
    end
end


N2 = next(N);

for j = J1:J2
    Nj = length(w{j});
    if radix2_flag == true
        t = (0:Nj-1)/Nj*N2/fs;
    else
        t = (0:Nj-1)/Nj*N/fs;
    end
    
    if 0
        Mj = 0.5/max(abs(w{j}));
    end
    
    %     if j <= J
    %         MM = Mj;
    %     else
    %         % j == J+1
    %         MM = 0.5 / max(abs(w{j}));
    %     end
    MM = Mj;    
    if j == J+1
        MM = min(Mj, 0.5 / max(abs(w{j})));
    end
    
    if STEM
        zz = zeros(size(t));
        if isreal(w{j})
            line([t; t],[zz; w{j}*MM]-j+J1-1, 'color','black' );
        else
            % line([t; t], [zz+w{j}*MM; zz+w{j}*MM]+J1-j-1, 'color','b')
            line([t; t],[zz; real(w{j})*MM]-j+J1-1,'color','black' );
            line([t; t],[zz; imag(w{j})*MM]-j+J1-1,'color','blue' );
        end
    else
        if isreal(w{j})
            line(t, w{j}*MM-j+J1-1, 'color', 'k');
        else
            line(t, real(w{j})*MM-j+J1-1, 'color','black');            
            line(t, imag(w{j})*MM-j+J1-1, 'color','blue');            
        end
    end
    
    if ENERGY
        % display fraction of total energy
        txt = sprintf('%5.2f%%',e(j)/sum(e)*100);
        th = text(1.02*N/fs,-j+J1-1,txt,'units','data');
        %         set(th, 'units', 'normalized')
    end
end
hold off
axis([0 N/fs -(J2-J1)-1.5 1])
title('SUBBANDS OF SIGNAL')

if fs == 1
    xlabel('TIME (SAMPLES)')
else
    xlabel('TIME (SECONDS)')
end
ylabel('SUBBAND')
set(gca,'Ytick',(J1-J2:0)-1,'YtickLabel',J2:-1:J1)

% Display parameter values in figure
if PARAM
    txt = sprintf('Q = %1.2f, r = %1.2f, Levels = %d',Q,r,J);
    % text(0,-0.05,txt,'units','normalized');

    tb = annotation('textbox',[.0 .0 .6 .05]);
    set(tb,...
        'string',txt,...
        'verticalalignment','bottom',...
        'edgecolor','none')

end

% txt = sprintf('Dilation = %.2f, Redundancy = %.2f',1/alpha, red);
% text(1,-0.05,txt,'units','normalized','horizontalalignment','right');

figure(gcf)

