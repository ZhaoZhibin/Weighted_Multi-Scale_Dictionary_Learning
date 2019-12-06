function [ Params ] = Config()
% This function performs parameter configuration
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6

%% Set the random seed to make sure the reproducibility
Params.random_seed = 25;          % The random state

%% Parameters of Generating Simulation
Params.Fs            = 10240;     % The sampling frequency of the simulation signal
Params.N             = 10240;      % The length of the signal
Params.mixture_ratio = [1.5, 1, 1]; % The mixing ratio of [impulses, harmonic, noise].
% Periodic Impulses
Params.AP         = 2;            % The amplitude of the impulse
Params.T          = 0.08;         % The fault period
Params.tau0       = 0.002;        % The initial phase
Params.f1         = 2050;         % The resonance frequency
Params.zeta       = 0.02;         % The left damping
% Harmonic interference
Params.Order = 2;                 % The Order of the interference
Params.CF    = 1000;               % The basic carrier frequency
Params.AM    = 180;                % The amplitude modulation frequency 
Params.FM    = 35;                % The frequency modulation frequency
Params.H     = [150, 550];        % The discrete frequencies
% noise type
Params.noise_type = 'Gaussian';   % The noise type can be 'Gaussian' or 'Laplacian'


%% Parameters of the HHLP
% TQWT parameters
Params.Q         = 8;            % The Q factor
Params.r         = 4;            % The redundant factor
Params.J         = 20;           % The level factor
% Dictionary learning parameters
Params.n         = 4;            % The size of each patch
Params.m         = 2*Params.n;   % The number of atoms
Params.init_E    = 26;           % The initialized threshold of OMP (This parameter needs to be fine-tuned)
Params.IterNum   = 3;            % Iterations of K-SVD
%Params.DictMS    = zeros(Params.n, Params.m);          % The initialized dictionary

end

