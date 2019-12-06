function [ Sig , t] = Generate_Simulation(Params)
% This function creates the simulation containing impulses, harmonic
% and noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params: a struct contains all the parameters
%%% Global parameters
%       Params.random_seed:   The random state
%       Params.Fs:            The sampling frequency of the simulation signal
%       Params.N:             The length of the signal
%       Params.mixture_ratio: The mixing ratio of [impulses, harmonic, noise]
%%% Periodic Impulses
%       Params.AP:            The amplitude of the impulse
%       Params.T:             The fault period
%       Params.f1:            The resonance frequency
%       Params.zeta:          The damping
%%% Harmonic interference
%       Params.Order:         The Order of the interference
%       Params.CF:            The basic carrier frequency
%       Params.AM:            The amplitude modulation frequency 
%       Params.FM:            The frequency modulation frequency
%       Params.H:             The discrete frequencies
%%% noise type
%       Params.noise_type:    The noise type can be 'Gaussian' or 'Laplacian'       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig:                        The generated signal
% t:                          The time index
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6
rng('default')
rng(Params.random_seed)
Fs = Params.Fs;
N = Params.N;
t = (0 : 1/Fs : (N-1)/Fs);
%% Generate the periodic impulses
AP = Params.AP;
f1 = Params.f1;
zeta = Params.zeta;          
T = Params.T;
NT = round(Fs * T);
Iter = floor(N / NT);
Sig_Impulse = [];
for i = 1 : Iter
    NT = round(Fs * (T + 0.001*(2*rand(1)-1)));
    t0 = (0 : 1/Fs : (NT-1)/Fs)';
    Sig_Impulse = [Sig_Impulse ; AP*exp(-zeta*2*pi*f1*t0).*sin(2*pi*f1*sqrt(1-zeta^2)*t0)];
end
if length(Sig_Impulse) >= N
    Sig_Impulse = Sig_Impulse(1:N);
else
    Sig_Impulse = [zeros(N-length(Sig_Impulse),1); Sig_Impulse];
end

%% Generate the Harmonic interference
Order = Params.Order;
CF = Params.CF;
AM = Params.AM;
FM = Params.FM;
H = Params.H;
Sig_Harmonic = 0;
for i = 1 : Order
    Sig_Harmonic = Sig_Harmonic + (1 + 0.5*cos(2*pi*i*AM*t)) .* cos(2*pi*i*CF*t + 0.5*cos(2*pi*i*FM*t));
end
for i = 1 : length(H)
    Sig_Harmonic = Sig_Harmonic + cos(2*pi*H(i)*t);
end

%% Generate the Noise
noise_type = Params.noise_type;
switch noise_type
    case 'Gaussian'
        Noise = randn(N , 1);
    case 'Laplacian'
        Noise = laprnd(N, 1);
    otherwise
        error('Unknown method.')
end
%% Generate the Combined Signal
mixture_ratio = Params.mixture_ratio;
Sig = mixture_ratio(1)*Sig_Impulse(:) + mixture_ratio(2)*Sig_Harmonic(:) + mixture_ratio(3)*Noise(:);
t = t(:);
end

