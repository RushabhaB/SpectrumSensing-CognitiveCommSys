function [CW_State, CW] = stage1_ED (nSU,CW_State,nSamples,E_s,N0,h_gain,th)

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

L = nSamples; % Number of sensing samples

for i = 1:length(CW_State)

for t = 1:length(th)
    mec = randsrc(nSU,L,[-1,1;0.5,0.5]);
    n =sqrt(N0./2)'.*(randn(nSU,L)+1i*randn(nSU,L)); % Gaussian noise, mean 0, variance 1
    h1 = (h_gain)*(randn(nSU,1)+1i*randn(nSU,1))./(sqrt(2)); % Generating Rayleigh channel coefficients
    h = repmat(h1,1,L); % Slow-fading
    if (CW_State(i) ==1)
        y = sqrt(E_s).*h.*mec + n; % Received signal y at the secondary user, abs(h) is the Rayleigh channel gain.
    else
        y = n;
    end
    energySU = mean(abs(y).^2,2);  % Received energy 
    CW(:,i) = double(energySU >= th(t)); % Checking whether the received energy is above the threshold
end
end
end
