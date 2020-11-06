% Cooperative Cognitive Radio Communication
% MATLAB code for assignment in AET G641 @ BITS Pilani
% Instructor: B. Sainath

% Students:
% Vandana Prasad - 2019H1240092P
% Rushabha Balaji - 2017A3PS0220P
% Vinay U Pai - 2017A3PS0131P

% Main program
clc;
clear all;

%Assumptions :
%Number of Secondary Users (SU) : 3
%SNR = 20dB for all signals from SU -> FC
%Rayleigh flat fading channel
% BPSK modulation

nSU = 10; % Number of Secondary Users
M = 2; %BPSK modulation
snr = 1:1:25;

for i = 1:1:25 %dB
%System Config
N0 = 1;  % Noise power
E_s = (10^(i/10))*N0; % Signal energy
%snr = 10*log10(E_s/N0); %SNR of the signal

% Rayleigh Fading Coefficients
H=(randn([nSU 1])+1j*randn([nSU 1]))/sqrt(2); % Channel Matrix

% Gaussian noise 
W=sqrt(N0)*(randn([nSU 1])+randn([nSU 1])*1j)/sqrt(2);% Noise vector of CSCG noise for SU
%W = [n0;n1;n2]; % Noise vector

% Data symbols

m = randi([0 1],[1,nSU])'; % Generating random message symbols
xb = sqrt(E_s)*exp(-1i*pi*2*(m)/M); % BPSK Data symbol for all the SU

X = diag(xb); % Generating the symbol matrix with the diagonal elements as the data symbols

Y = X*H + W; % Output symbols at the Fusion Centre (FC)

%Channel  Estimation
H_LS = inv(X*X')*X'*Y; %MISO Least Square detection using pilots

% Detection using estimated channel 
bpsk_demod = comm.BPSKDemodulator; %Creating the bpsk demodulator object

% Detecting the received bits from the ZF equalisation of received vector
xb_det = inv(diag(H_LS))*Y;
m_det = bpsk_demod(xb_det);

%BER
[nErr(floor(i/1)),ratio(floor(i/1))] = biterr(m,m_det);
end

plot(snr,nErr)