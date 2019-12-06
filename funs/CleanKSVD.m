function [y_out] = CleanKSVD(y, Params)
% This function realizes signal denoising using KSVD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y:                          The input signal
% Params: a struct contains all the parameters
%%% Dictionary learning parameters
%       Params.n:             The size of each patch
%       Params.m:             The number of atoms
%       Params.E:             The threshold of OMP
%       Params.IterNum:       Iterations of K-SVD
%       Params.DictMS:        The initialized dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_out:                      The denoised signal


% Reference: 'A weighted multi-scale dictionary learning model and its applications on bearing fault diagnosis', Journal of Sound and Vibration, 2019
% HomePages: https://zhaozhibin.github.io/
% Author   : Zhibin Zhao
% Place    : Xi'an Jiaotong University
% Email    : zhibinzhao1993@gmail.com
% Date     : 2017.2
if nargin < 2
    Params = struct([]);
end
% initialization
if ~isfield(Params, 'n')
    n = 200;
else
    n = Params.n;
end
if ~isfield(Params, 'm')
    m = n*2;
else
    m = Params.m;
end
if ~isfield(Params, 'E')
    E = 10;
else
    E = Params.E;
end
if ~isfield(Params, 'IterNum')
    IterNum = 20;
else
    IterNum = Params.IterNum;
end
if ~isfield(Params, 'DictMS')
    initdict = odctdict(n , m);
else
    initdict = Params.DictMS;
end

y = y(:);
N = length(y);
Olop = 1;                                                            % the overlapping coefficient
% Extract the patches
NSample = floor((N - n)/Olop)+1;                     % the number of the patches
y = y(1 : (NSample-1)*Olop+n);
Y = zeros(n , NSample);
% Put the patches into the predefined matrix
for i = 1 : NSample
    Start  = 1+(i-1)*Olop;
    patch  = y(Start : Start+n-1);
    Y(:,i) = patch(:);   
end

params.Edata = E;
params.iternum = IterNum;
params.data = Y;
params.initdict = initdict;

[Dict, X] = ksvd(params, '');                               % Train the dictionary

Yd = Dict*X;

% Solve the formula : argmin sum(||Dict*X - R*y||^2)
Weight = zeros(length(y),1);
y_out     = zeros(length(y),1);
for i = 1 : NSample
    Start  = 1+(i-1)*Olop;
    patch  = Yd(:,i);
    y_out(Start : Start+n-1) = y_out(Start : Start+n-1) + patch;
    Weight(Start : Start+n-1) = Weight(Start : Start+n-1) + 1;
end
y_out = y_out ./ Weight;
end

