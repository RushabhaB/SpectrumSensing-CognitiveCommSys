function [CW_State, CW] = stage1_ED (nSU,nCodeword,nSamples,E_s,fa)
CW=[];

CW_State = randi([0,1], [1,nCodeword]);

L = nSamples; % Number of sensing samples
iter =10^5; % Number of iterations 
%fa =0.01:0.03:1; % Probability of False Alarm
%fa = 0.05;
E_p_l = E_s; % SNR in linear scale
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
thresh_sim(t) = energy_desc(ceil(fa(t)*iter)); % Threshold obtained by simulations
th(t) = (qfuncinv(fa(t)/2))^2*(N0/2);
end

for i = 1:length(CW_State)

for t = 1:length(fa)
    s = [];  
    h = [];
    mes = randi([0 1],nSU,L); % Generating 0 and 1 with equal probability for BPSK
    n =sqrt(N0)*randn(nSU,L); % Gaussian noise, mean 0, variance 1
    s = (2.*(mes)-1); % BPSK modulation
    h1 = (randn(nSU,1)+1i*randn(nSU,1))./(sqrt(2)); % Generating Rayleigh channel coefficients
    h = repmat(h1,1,L); % Slow-fading
    if (CW_State(i) ==1)
        y = sqrt(E_p_l).*abs(h).*s + n; % Received signal y at the secondary user, abs(h) is the Rayleigh channel gain.
    else
        y = n;
    end
    energy = (abs(y).^2);  % Received energy 
    
for k=1:nSU % Number of Monte Carlo Simulations
    energy_fin_1(k,:) =sum(energy(k,:)); % Test Statistic for the energy detection 
end

    detect = ((energy_fin_1/10) >= th(t)); % Checking whether the received energy is above the threshold
    %i_prime = sum(detect); % Count how many times out of 'iter', the received energy is above the threshold.
    CW(i,:)=double(detect);
    
    %final(final==0)=-1; % resultant matrix transmitted by SU
end
end
end
