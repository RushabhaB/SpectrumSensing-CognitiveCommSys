function [CW_State, CW] = stage1_ED (nSU,nCodeword,nSamples,E_s,N0,h_gain,th)

% This function provides the results of the status of the PU from all the
% SUs based on the Energy Detection (ED) technique.
%  Inputs to the function :
%        nSU = number of Users
%        nCodeWords = number of codewords
%        nSamples = number of samples for which the channel remains
%        constant
%        E_s = energy of the symbol
%        N0 = noise power
%        h_gain = Channel power gain
%        fa = local false alarm rate vector
% Outputs the function provides :
%        CW_State : Actual state of the PU whose size is 1xnCodeWords
%        CW : Detected PU state by each SU , matrix size is nCodeWordsxnSU
 
CW=[];

CW_State = randi([0,1], [1,nCodeword]);

L = nSamples; % Number of sensing samples
iter =10^5; % Number of iterations 

E_p_l = E_s; % SNR in linear scale

%th = (qfuncinv(fa/2))^2*(N0/2);

for i = 1:length(CW_State)

for t = 1:length(th)
    s = [];  
    h = [];
    mes = randi([0 1],nSU,L); % Generating 0 and 1 with equal probability for BPSK
    n =sqrt(N0').*randn(nSU,L); % Gaussian noise, mean 0, variance 1
    s = (2.*(mes)-1); % BPSK modulation
    h1 = (h_gain)*(randn(nSU,1)+1i*randn(nSU,1))./(sqrt(2)); % Generating Rayleigh channel coefficients
    h = repmat(h1,1,L); % Slow-fading
    if (CW_State(i) ==1)
        y = sqrt(E_p_l)'.*abs(h).*s + n; % Received signal y at the secondary user, abs(h) is the Rayleigh channel gain.
    else
        y = n;
    end
    energy = (abs(y).^2);  % Received energy 
    
for k=1:nSU % Number of Monte Carlo Simulations
    energy_fin_1(k,:) =sum(energy(k,:)); % Test Statistic for the energy detection 
end

    detect = ((energy_fin_1/L) >= th(t)); % Checking whether the received energy is above the threshold
    %i_prime = sum(detect); % Count how many times out of 'iter', the received energy is above the threshold.
    CW(i,:)=double(detect);
    
    %final(final==0)=-1; % resultant matrix transmitted by SU
end
end
end
