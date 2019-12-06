%% Demo1: 
% HomePages: https://zhaozhibin.github.io/
% Author   : Zhibin Zhao
% Place    : Xi'an Jiaotong University
% Email    : zhibinzhao1993@gmail.com
% Date     : 2019.12
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));

Params = Config();
rng('default')
rng(Params.random_seed)
[ Sig , t] = Generate_Simulation(Params);


%% Estimate the energy of the noise
NoiseSigma = NoiseEstimate(Sig);
%% Perform weighted multi-scale dictionary learning
Params.init_E = Params.init_E * NoiseSigma;
tic
[y_WMSDL] = WMSDL(Sig, Params);
time_WMSDL = toc;
%% Perform KSVD denoising
Params.n = 200;
Params.m = Params.n * 2;
Params.E = 19.5 * NoiseSigma * getConstant(Params.n);
tic
[y_KSVD] = CleanKSVD(Sig, Params);
time_KSVD = toc;
%% Perform the square envelope spectrum (SES)
[ yf_WMSDL, ~ ] = Hilbert_envelope( y_WMSDL , Params.Fs , 1);
[ yf_KSVD, f ] = Hilbert_envelope( y_KSVD , Params.Fs , 1);
%% Plot the results
figure(1)
subplot(311)
plot(t, Sig)
title('Original Signal')
ylabel('Amp (g)')
subplot(312)
plot(t, y_WMSDL)
title(['WMSDL, Computational time = ' num2str(time_WMSDL) 's'])
ylabel('Amp (g)')
subplot(313)
plot(t, y_KSVD)
title(['KSVD, Computational time = ' num2str(time_KSVD) 's'])
ylabel('Amp (g)')
xlabel('Time (s)')
filename = ['results', filesep, sprintf('Demo1_Extracted_Time.pdf')];
print(filename, '-dpdf');
figure(2)
subplot(211)
plot(f(1:800), yf_WMSDL(1:800))
title(['WMSDL, Computational time = ' num2str(time_WMSDL) 's'])
ylabel('Amp (g)')
subplot(212)
plot(f(1:800), yf_KSVD(1:800))
title(['KSVD, Computational time = ' num2str(time_KSVD) 's'])
ylabel('Amp (g)')
xlabel('Frequency (Hz)')
filename = ['results', filesep, sprintf('Demo1_Extracted_SES.pdf')];
print(filename, '-dpdf');
% If you feel our HHLP is useful for your research, please consider citing our paper:
% @article{zhao2019weighted,
%   title={A weighted multi-scale dictionary learning model and its applications on bearing fault diagnosis},
%   author={Zhao, Zhibin and Qiao, Baijie and Wang, Shibin and Shen, Zhixian and Chen, Xuefeng},
%   journal={Journal of Sound and Vibration},
%   volume={446},
%   pages={429--452},
%   year={2019},
%   publisher={Elsevier}
% }