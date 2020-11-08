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
bpsk_demod = comm.BPSKDemodulator; %Creating the bpsk demodulator object

nSU = 3; % Number of Secondary Users
M = 2; %BPSK modulation
%snr = 1:1:25; % Use this array to generate ber vs snr plot

%for i = 1:1:25 %dB

%System Config
N0 = 1;  % Noise power
nCodeWords = 500;
E_s =100;
%E_s = (10^(i/10))*N0; % Signal energy
snr = 10*log10(E_s/N0); %SNR of the signal
nSamples = 3 ; % Number of samples across which we assume the H_coeff to be constant
m_p = ones([nSU 1]); % Pilot bits
CW = zeros([nSU,nCodeWords]); % Actual Codewords received
CW_det = zeros([nSU,nCodeWords]);% Detected codewords
% Rayleigh Fading Coefficients
H=(randn([nSU nCodeWords])+1j*randn([nSU nCodeWords]))/sqrt(2); % Channel Matrix

% Gaussian noise 
W=sqrt(N0)*(randn([nSU nCodeWords])+randn([nSU nCodeWords])*1j)/sqrt(2);% Noise vector of CSCG noise for SU

for i = 1: nSamples+1 :nCodeWords
CW(:,i) = m_p;
CW_det(:,i) = m_p;
x_p = sqrt(E_s)*exp(-1i*pi*2*(m_p)/M); % BPSK Pilot symbol for all the SU

X_p = diag(x_p); % Generating the symbol matrix with the diagonal elements as the data symbols

Y_p = X_p*H(:,i) + W(:,i); % Output symbols at the Fusion Centre (FC)

%Channel  Estimation
H_LS_p = inv(X_p*X_p')*X_p'*Y_p; %MISO Least Square detection using pilots

for j = i+1:i+nSamples

% Data symbols
m = randi([0 1],[1,nSU])'; % Generating random message symbols
xb = sqrt(E_s)*exp(-1i*pi*2*(m)/M); % BPSK Data symbol for all the SU
X_b = diag(xb); % Generating the symbol matrix with the diagonal elements as the data symbols
Y_b = X_b*H(:,j) + W(:,j); % Output symbols at the Fusion Centre (FC)

% Detection using estimated channel 

% Detecting the received bits from the ZF equalisation of received vector
xb_det = inv(diag(H_LS_p))*Y_b;
m_det = bpsk_demod(xb_det);

%BER
%[nErr(floor(i/1)),ratio(floor(i/1))] = biterr(m,m_det); %Generating the
%number of errors array and ratio array ( include in snr for loop)
[nErr(j),ratio(j)] = biterr(m,m_det); %Generating the number of errors array and ratio array (without for loop)

CW(:,j) = m;
CW_det(:,j)= m_det;
end
end

[p_md,p_fa] = md_fa(CW,CW_det,nSamples,nCodeWords);

p_md
p_fa

