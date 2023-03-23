%% A comparison of OFDM Radar signal processing methods
% Description: This project compares some classical sensing algorithm.
% Author: Yunbo HU (SIMIT, UCAS)
clc;
clear;
addpath('fig');
global c0 fc lambda M N delta_f Ts
%% ISAC Transmitter
% System parameters
c0 = 3e+8;  % velocity of light
fc = 30e+9; % carrier frequency
lambda = c0 / fc; % wavelength
M = 1024; % number of subcarriers
N = 15; % number of symbols per subframe
delta_f = 120e+3; % subcarrier spacing
T = 1 / delta_f; % symbol duration
Tcp = T / 4; % cyclic prefix duration
Ts = T + Tcp; % total symbol duration
CPsize = M / 4; % cyclic prefix length
bitsPerSymbol = 4; % bits per symbol
qam = 2^(bitsPerSymbol); % 16-QAM modulation

% Transmit data
data = randi([0 qam - 1], M, N);
TxData = qammod(data, qam, 'gray');

% OFDM modulator
TxSignal = ifft(TxData, M); % IFFT
TxSignal_cp = [TxSignal(M - CPsize + 1: M, :); TxSignal]; % add CP
TxSignal_cp = reshape(TxSignal_cp, [], 1); % time-domain transmit signal

%% Channel
% Sensing Data Generation
SNR = 10;
RxSignal = sensingSignalGen(TxSignal_cp, 30,20,SNR);


%% OFDM Radar Receivers

% 1. 2DFFT based   (classical method, reference is omitted here)

Rx = RxSignal(1:size(TxSignal_cp,1),:); 

Rx = reshape(Rx,[],N);
Rx = Rx(CPsize + 1 : M + CPsize,:); % remove CP and reshape it into a matrix
Rx_dem = fft(Rx,M);
CIM_2dfft = Rx_dem .* conj(TxData); %elementwise-multiply(equals to match filtering)

RDM_2dfft = fft(ifft(CIM_2dfft,M).',10*N);

% plot the range doppler map
figure(1);
range_2dfft = linspace(0,c0/(2*delta_f),M+1);
range_2dfft = range_2dfft(1:M);

velocity_2dfft = linspace(0,lambda/2/Ts,10*N+1);
velocity_2dfft = velocity_2dfft(1:10*N);

[X,Y] = meshgrid(range_2dfft,velocity_2dfft);
RDM_2dfft_norm = 10*log10( abs(RDM_2dfft) / max(abs(RDM_2dfft),[],'all'));
surf(X,Y,(RDM_2dfft_norm));
title('2D-FFT based method');
xlabel('range(m)');
ylabel('velocity(m/s)');
savefig('fig/figure1.fig');
% 2. CCC-based     (Method proposed by Kai Wu et al.)
figure(2);
% setting parameters for CCC-based sensing method
% please refer to the paper for the meaning of these paramaters
mildM = 512;
Qbar = 64;
mildQ = 128;
% CCC
[r_cc,RDM] = cccSensing(RxSignal,TxSignal_cp,mildM,Qbar,mildQ);   % ccc sensing

%plot the range doppler map
Tsa = 1/delta_f/M;
mildN = floor((length(TxSignal_cp)-Qbar-mildQ)/(mildM - Qbar));
range_ccc = linspace(0,c0/2*Tsa*mildM, mildM+1);
doppler_ccc = linspace(0,lambda/(mildM-Qbar)/Tsa/2,10*mildN+1);
range_ccc = range_ccc(1:mildM);
doppler_ccc = doppler_ccc(1:10*mildN);

RDM_norm = 10*log10(abs(RDM)/max(abs(RDM),[],'all'));
[X,Y] = meshgrid(range_ccc,doppler_ccc);
surf(X,Y,(RDM_norm)); % plot the range-doppler map
title('CCC based method');
xlabel('range(m)');
ylabel('velocity(m/s)');
savefig('fig/figure2.fig');

% 3. Super resolution sensing method
% 3.1 MUSIC based (a time consuming but precise method)
CIM = Rx_dem ./(TxData); 

[P_music_range,P_music_velo] = MUSICforOFDMsensing(CIM,1);


% plot the MUSIC power spectrum
figure(3);
title('MUSIC for OFDM sensing');
subplot(1,2,1);
plot(linspace(0,100,M),abs(P_music_range)/max(abs(P_music_range)));
ylabel('Pmusic');
xlabel('range(m)');
ylim([10^-3,1]);
title('MUSIC for range estimation');
subplot(1,2,2);
plot(linspace(0,100,M),abs(P_music_velo)/max(abs(P_music_velo)));
ylabel('Pmusic');
xlabel('velocity(m/s)');
ylim([10^-3,1]);
title('MUSIC for velocity estimation');

savefig('fig/figure3.fig');

% 3.2 ESPRIT based method
[range,velocity] = ESPRITforOFDMsensing(CIM,1);
fprintf('The estimation result of TLS-ESPRIT is :\n');
fprintf('Range = %f\n',range);
fprintf('Velocity = %f\n',velocity);

