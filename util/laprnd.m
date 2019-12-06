function y  = laprnd(m, n, mu, sigma)
% This function generates i.i.d. laplacian noise
%   [m, n]  : the dimension of y.
%   mu      : The mean of the noise (Default: 0)
%   sigma   : The standard deviation of the noise (Default: 1)

% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6


if nargin < 4
    sigma = 1;
end
if nargin < 3
    mu = 0;
end
if nargin < 2
    error('You must input the dimension of the generated noise')
end
% Generate Laplacian noise
a = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(a).* log(1- 2* abs(a));