clc
close all
clear all
L = 10; % Number of sensing samples
iter =10^5; % Number of iterations 
%fa =0.01:0.03:1; % Probability of False Alarm
fa = 0.01;
snr_db = 0; % SNR in dB
snr = 10.^(snr_db./10); % SNR in linear scale
for t = 1:length(fa) % Calculating threshold for each value of false alarm
 energy_fin = []; 
 n = [];
 y = [];
 energy = [];
 n=randn(iter,L); % Gaussian noise, mean 0, variance 1
 y = n; % Received signal at the secondary user
 energy = abs(y).^2; % Energy of received signal over L sensing samples 
for k=1:iter 
 energy_fin(k,:) =sum(energy(k,:)); % Test Statistic of the energy detection
end
energy_fin = energy_fin'; % Taking transpose to arrage values in descending order
energy_desc = sort(energy_fin,'descend')'; % Arrange values in descending order
thresh_sim(t) = energy_desc(ceil(fa(t)*iter)); % Threshold obtained by simulations
end

for t = 1:length(fa)

    s = [];  
    h = [];
    mes = randi([0 1],iter,L); % Generating 0 and 1 with equal probability for BPSK
    s = (2.*(mes)-1); % BPSK modulation
    h1 = (randn(iter,1)+j*randn(iter,1))./(sqrt(2)); % Generating Rayleigh channel coefficients
    h = repmat(h1,1,L); % Slow-fading 
    y = sqrt(snr).*abs(h).*s + n; % Received signal y at the secondary user, abs(h) is the Rayleigh channel gain.
    energy = (abs(y).^2);  % Received energy 
    
for k=1:iter % Number of Monte Carlo Simulations
    energy_fin(k) =sum(energy(k,:)); % Test Statistic for the energy detection 
end

    detect = (energy_fin >= thresh_sim(t)); % Checking whether the received energy is above the threshold
    i_prime = sum(detect); % Count how many times out of 'iter', the received energy is above the threshold.
    final=double(detect);
    final(final==0)=-1; % resultant matrix transmitted by SU
end
i_prime