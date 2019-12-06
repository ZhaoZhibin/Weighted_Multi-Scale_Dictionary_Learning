function [ y_out ] = WMSDL( y , Params)
% This function solves the weighted multi-scale dictionary learning problem:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y:                          The input signal
% Params: a struct contains all the parameters
%%% TQWT parameters
%       Params.Q:             The Q factor
%       Params.r:             The redundant factor (default: 5)
%       Params.J:             The level factor (default: 10)
%%% Dictionary learning parameters
%       Params.n:             The size of each patch
%       Params.m:             The number of atoms
%       Params.init_E:        The initialized threshold of OMP (This parameter needs to be fine-tuned)
%       Params.IterNum:       Iterations of K-SVD
%       Params.DictMS:        The initialized dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_out:                      The denoised signal


% Reference: 'A weighted multi-scale dictionary learning model and its applications on bearing fault diagnosis', Journal of Sound and Vibration, 2019
% Homepage : https://zhaozhibin.github.io/
% Author   : Zhibin Zhao
% Place    : Xi'an Jiaotong University
% Email    : zhibinzhao1993@gmail.com
% Date     : 2017.2

if nargin < 2
    Params = struct([]);
end
% initialization
if ~isfield(Params, 'Q')
    Q = 8;
else
    Q = Params.Q;
end
if ~isfield(Params, 'r')
    r = 4;
else
    r = Params.r;
end
if ~isfield(Params, 'J')
    J = 20;
else
    J = Params.J;
end
if ~isfield(Params, 'init_E')
    init_E = 10;
else
    init_E = Params.init_E;
end

N = length(y);
Params.N = N;
Norm = ComputeNow(N, Q, r, J);
w = tqwt(y, Q, r, J);
d = tqwt(zeros(size(y)),Q,r,J);
for S = J+1 : -1 : 1   
    coefd = w{S};
    % Calculate the weight
    d_tmp = d;
    d_tmp{S} = coefd;
    y_single = itqwt(d_tmp, Q, r, N);
    Ks = kurtosis(y_single);
    Params.E = init_E * Norm(S) / Ks;
    % Dictionary learning of approximate coefficients 
    [dcoefd] = CleanKSVD(coefd, Params);
    w{S} = dcoefd';  
end
y_out = itqwt(w, Q, r, N);
end


