% Cooperative Cognitive Radio Communication
% MATLAB code for assignment in AET G641 @ BITS Pilani
% Instructor: B. Sainath

% Students:
% Vandana Prasad - 2019H1240092P
% Rushabha Balaji - 2017A3PS0220P
% Vinay U Pai - 2017A3PS0131P

% Function for Stage 1
% Refer to stage1_ED.m for functionality
% Only change is that E_s is being varied while Pfa is constant

function [CW_State, CW] = stage1_ED_Es_change (nSU,nCodeword,nSamples,E_s,fa)

%CW_State is the ground truth (1s mean active, 0s mean idle)
%CW contains the codewords received by the FC from the SUs

CW=[];

CW_State = randi([0,1], [1,nCodeword]); % Ground truth on whether channel is active (1) or idle (0) randomly

L = nSamples; % Coherence interval
iter =10^5; % Number of iterations 

for t = 1:length(fa) % Calculating threshold for each value of false alarm
 energy_fin = []; 
 N0=1; % Noise power
 n = [];
 y = [];
 energy = [];
 n=sqrt(N0)*randn(iter,L); % Gaussian noise, mean 0, variance 1
 y = n; % Received signal at the secondary user
 energy = abs(y).^2; % Energy of received signal over L sensing samples 
for k=1:iter 
 energy_fin(k,:) = sum(energy(k,:)); % Test Statistic of the energy detection
end
energy_fin = energy_fin'; % Taking transpose to arrage values in descending order
energy_desc = sort(energy_fin,'descend')'; % Arrange values in descending order
thresh_sim = energy_desc(ceil(fa*iter)); % Threshold obtained by simulations
th = (qfuncinv(fa/2))^2*(N0/2); %Getting threshold from required local P_fa
end

for i = 1:length(CW_State)

for t = 1:length(E_s)
    s = [];  
    h = [];
    mes = randi([0 1],nSU,L); % Generating 0 and 1 with equal probability for BPSK
    n =sqrt(N0)*randn(nSU,L); % Gaussian noise, mean 0, variance N0
    s = (2.*(mes)-1); % BPSK modulation
    h1 = (randn(nSU,1)+1i*randn(nSU,1))./(sqrt(2)); % Generating Rayleigh channel coefficients
    h = repmat(h1,1,L); % Channel is constant over coherence interval
    if (CW_State(i) ==1)
        y = sqrt(E_s(t)).*abs(h).*s + n; % Received signal y at the SU if channel is active, abs(h) is the Rayleigh channel gain.
    else
        y = n;  % Received signal y is only noise isf channel is idle
    end
    energy = (abs(y).^2);  % Received energy 
    
for k=1:nSU % Number of SUs
    energy_fin_1(k,:) =sum(energy(k,:)); % Test Statistic for the energy detection 
end

    detect = ((energy_fin_1/L) >= th); % Checking whether the received energy is above the threshold
    CW(i,:)=double(detect); % Convert all above threshold to 1 and below to 0
    

end
end
end
