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

%Reference: "Cooperative Spectrum Sensing Using Maximum a Posteriori 
%as a Detection Technique for Dynamic Spectrum Access Networks,"

%Assumptions :
%Number of Secondary Users (SU) : 3
%Rayleigh flat fading channel
%BPSK modulation
bpsk_demod = comm.BPSKDemodulator; %Creating the bpsk demodulator object

nSU = 3; % Number of Secondary Users (Keep this as 3 to get meaning ful MAP results)
M = 2; %BPSK modulation
h_gain = 1;

%System Config
N0 = 10;  % Noise power
nCodeWords = 10^4; % Number of codewords per iteration
nSamples = 4 ; % Number of samples across which we assume the H_coeff to be constant
m_p = ones([nSU 1]); % Pilot bits

power_noise_ratio_db= 0:2:50; % Range of SNR (in dB) over which we want to plot the parameters
power_noise_ratio = 10.^(power_noise_ratio_db/10); % In linear scale
E_s = power_noise_ratio.*N0; % Signal power

fa = 0.05; % Desired local SU P_fa
[th,CW_p] = MAP_est_Es_change(fa,N0,E_s,nSU); % Get threshold and probabliity map for each codeword


iter = 100; %Number of MC simulations

% Initialising the modulated BPSK symbols for all codewords
allCW = [1,1,1;1,1,0;1,0,1;1,0,0;0,1,1;0,1,0;0,0,1;0,0,0]'; 


for r=1:iter
    r
%Majority arrays
p_md_arr = [];
p_fa_arr = [];
p_md_mmse_arr = [];
p_fa_mmse_arr = [];
p_md_ideal_arr=[];
p_fa_ideal_arr =[];

% MAP arrays
p_md_MAP_ideal_arr = [];
p_fa_MAP_ideal_arr = [];
p_md_MAP_LS_arr = [];
p_fa_MAP_LS_arr = [];
p_md_MAP_mmse_arr = [];
p_fa_MAP_mmse_arr = [];

for k=1:length(E_s)
[CW,CW_detSU] = stage1_ED_Es_change (nSU,nCodeWords,nSamples,E_s(k),fa); %Get ground truth for channel status and codewords at FC
x_allCW = sqrt(E_s(k)).*exp(-1i*pi*2*(allCW)/M);
%Initialise arrays for Ideal (Perfect CSI, MMSE and LS)
CW_ideal = ones(size(CW));
CW_LS = ones(size(CW));
CW_mmse = ones(size(CW));

pilot_loc = 1:nSamples+1:nCodeWords ; % Location of pilot bits (every 5th bit as coherence interval is 5 samples)
for i = 1:nSamples+1:nCodeWords
    CW(i) = m_p(1);%The pilots are inserted here as well, though of no consequence it just makes it easier to compare and avoids confusion
    CW_detSU(:,i) = m_p;  %Inserting the pilot vectors into codeword by forcing  the every nSample'th bits to be the pilot
end
 
% Rayleigh Fading Coefficients
H=(randn([nSU ceil(nCodeWords/(nSamples+1))])+1j*randn([nSU ceil(nCodeWords/(nSamples+1))]))/sqrt(2); % Channel Matrix

% Gaussian noise 
W=sqrt(N0)*(randn([nSU nCodeWords])+randn([nSU nCodeWords])*1j)/sqrt(2);% Noise vector of CSCG noise for SU

%Data symbols
for i = 1: nSamples+1 :nCodeWords
x_p = sqrt(E_s(k)).*exp(-1i*pi*2*(m_p)/M); % BPSK Pilot symbol for all the SU

X_p = diag(x_p); % Generating the symbol matrix with the diagonal elements as the data symbols

Y_p = X_p*H(:,ceil(i/(nSamples+1))) + W(:,i); % Output symbols at the Fusion Centre (FC)

%Channel  Estimation
H_LS_p(:,ceil(i/(nSamples+1))) = inv(X_p*X_p')*X_p'*Y_p; %MISO Least Square detection using pilots
u = (1./(E_s(k) + N0)).*x_p; %E[h^2] is  1. Also ||x||^2 is E_s
H_mmse_p(:,ceil(i/(nSamples+1))) = conj(u).*Y_p; %Standard formula from Tse Vishwanath

end



% Determine the data 
for i = 1:nSamples+1:nCodeWords
for j = i+1:i+nSamples

% Data symbols
m = CW_detSU(:,j); % Transmitting the received PU signal by the SU to the FC
xb = sqrt(E_s(k)).*exp(-1i*pi*2*(m)/M); % BPSK Data symbol for all the SU
X_b = diag(xb); % Generating the symbol matrix with the diagonal elements as the data symbols
Y_b = X_b*H(:,ceil(i/(nSamples+1))) + W(:,j); % Output symbols at the Fusion Centre (FC)

% Detection using estimated channel 

% Detecting the received bits from the ZF equalisation (Perfect CSI) of received vector
xb_det_ideal = inv(diag(H(:,ceil(i/(nSamples+1)))))*Y_b;
m_det_ideal =  bpsk_demod(xb_det_ideal);

% Detecting the received bits from the ZF equalisation (LS) of received vector
xb_det = inv(diag(H_LS_p(:,ceil(i/(nSamples+1)))))*Y_b;
m_det = bpsk_demod(xb_det);

% Detecting the received bits from the ZF equalisation (MMSE) of received vector
xb_det_mmse = inv(diag(H_mmse_p(:,ceil(i/(nSamples+1)))))*Y_b;
m_det_mmse = bpsk_demod(xb_det_mmse);

%Begin the MAP estimate 

%Reference Algorithm 1

Ka_min=zeros(1,length(allCW));
Ki_min=zeros(1,length(allCW));
Ka_min_LS=zeros(1,length(allCW));
Ki_min_LS=zeros(1,length(allCW));
Ki_min_mmse=zeros(1,length(allCW));
Ka_min_mmse=zeros(1,length(allCW));
for q=1:length(allCW)
    % Ideal Channel combiner
    Ka_min(q) = sum(abs(Y_b-diag(x_allCW(:,q))*H(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(1,q,k));
    Ki_min(q)= sum(abs(Y_b-diag(x_allCW(:,q))*H(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(2,q,k));

    % LS channel combiner
    Ka_min_LS(q) = sum(abs(Y_b-diag(x_allCW(:,q))*H_LS_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(1,q,k));
    Ki_min_LS(q) = sum(abs(Y_b-diag(x_allCW(:,q))*H_LS_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(2,q,k));

    % MMSE channel combiner
    Ka_min_mmse(q) = sum(abs(Y_b-diag(x_allCW(:,q))*H_mmse_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(1,q,k));
    Ki_min_mmse(q) = sum(abs(Y_b-diag(x_allCW(:,q))*H_mmse_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(2,q,k));
end

if (min(Ka_min) < min(Ki_min))
    CW_ideal(j) = 1;
else
    CW_ideal(j) = 0;
end
    
if (min(Ka_min_LS) < min(Ki_min_LS))
    CW_LS(j) = 1;
else
    CW_LS(j) = 0;
end

if (min(Ka_min_mmse) < min(Ki_min_mmse))
    CW_mmse(j) = 1;
else
    CW_mmse(j) = 0;
end

%BER

[nErr(j),ratio(j)] = biterr(m,m_det); %Generating the number of errors array and ratio array for ls (without for loop)
[nErr_mmse(j),ratio_mmse(j)] = biterr(m,m_det_mmse); %Generating the number of errors array and ratio array for mmse (without for loop)

CW_detFC_ideal(:,j) = m_det_ideal; %Detected codewords (perfect CSI) after ZF 
CW_detFC(:,j)= m_det;  %Detected codewords (LS) after ZF 
CW_detFC_mmse(:,j)= m_det_mmse;  %Detected codewords (MMSE) after ZF 
end
end

%Get Pmd and Pfa for Majority combining
[p_md_ideal,p_fa_ideal] = md_fa(CW,CW_detFC_ideal,nSamples,nCodeWords);
[p_md,p_fa] = md_fa(CW,CW_detFC,nSamples,nCodeWords);
[p_md_mmse,p_fa_mmse] = md_fa(CW,CW_detFC_mmse,nSamples,nCodeWords);

%Get Pmd and Pfa for MAP combining

[p_md_MAP_ideal,p_fa_MAP_ideal] = md_fa_MAP(CW,CW_ideal,nSamples,nCodeWords);
[p_md_MAP_LS,p_fa_MAP_LS] = md_fa_MAP(CW,CW_LS,nSamples,nCodeWords);

[p_md_MAP_mmse,p_fa_MAP_mmse] = md_fa_MAP(CW,CW_mmse,nSamples,nCodeWords);

% end MAP estimation

% Store the Majority combining values in an array for every value of P/N0
p_md_ideal_arr = [p_md_ideal_arr p_md_ideal];
p_fa_ideal_arr = [p_fa_ideal_arr p_fa_ideal];

p_md_arr = [p_md_arr p_md];
p_fa_arr = [p_fa_arr p_fa];

p_md_mmse_arr = [p_md_mmse_arr p_md_mmse];
p_fa_mmse_arr = [p_fa_mmse_arr p_fa_mmse];


% Store the MAP combining values in an array for every value of P/N0
p_md_MAP_ideal_arr = [p_md_MAP_ideal_arr p_md_MAP_ideal];
p_fa_MAP_ideal_arr = [p_fa_MAP_ideal_arr p_fa_MAP_ideal];
% 
p_md_MAP_LS_arr = [p_md_MAP_LS_arr p_md_MAP_LS];
p_fa_MAP_LS_arr = [p_fa_MAP_LS_arr p_fa_MAP_LS];
% 
p_md_MAP_mmse_arr = [p_md_MAP_mmse_arr p_md_MAP_mmse];
p_fa_MAP_mmse_arr = [p_fa_MAP_mmse_arr p_fa_MAP_mmse];


end

%Values of the parameters stored for each iteration
p_md_ideal_array(r,:) = p_md_ideal_arr;
p_fa_ideal_array(r,:) = p_fa_ideal_arr;

p_md_LS_array(r,:) = p_md_arr;
p_fa_LS_array(r,:) = p_fa_arr;

p_md_mmse_array(r,:) = p_md_mmse_arr;
p_fa_mmse_array(r,:) = p_fa_mmse_arr;

p_md_MAP_ideal_array(r,:) = p_md_MAP_ideal_arr;
p_fa_MAP_ideal_array(r,:) = p_fa_MAP_ideal_arr;
p_md_MAP_LS_array(r,:) = p_md_MAP_LS_arr;
p_fa_MAP_LS_array(r,:) = p_fa_MAP_LS_arr;
p_md_MAP_mmse_array(r,:) = p_md_MAP_mmse_arr;
p_fa_MAP_mmse_array(r,:) = p_fa_MAP_mmse_arr;

end

%Get mean across number of iterations
p_md_ideal_final = mean(p_md_ideal_array,1);
p_fa_ideal_final = mean(p_fa_ideal_array,1);

p_md_LS_final = mean(p_md_LS_array,1);
p_fa_LS_final = mean(p_fa_LS_array,1);

p_md_mmse_final = mean(p_md_LS_array,1);
p_fa_mmse_final = mean(p_fa_LS_array,1);

p_md_MAP_ideal_final = mean(p_md_MAP_ideal_array,1);
p_fa_MAP_ideal_final = mean(p_fa_MAP_ideal_array,1);

p_md_MAP_LS_final = mean(p_md_MAP_LS_array,1);
p_fa_MAP_LS_final = mean(p_fa_MAP_LS_array,1);

p_md_MAP_mmse_final = mean(p_md_MAP_mmse_array,1);
p_fa_MAP_mmse_final = mean(p_fa_MAP_mmse_array,1);


%Plot the figures 

% MAP vs Majority Ideal
 figure(1)
 grid on 
 plot(power_noise_ratio_db,p_md_ideal_final ,'k-d','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,p_md_MAP_ideal_final ,'r-d','LineWidth',2);
 hold on
 xlabel('$P/N_0$ (dB)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of misdetection ($P_{md}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('Majority Combiner', 'MAP Combiner','Location','northeast','FontSize',10,'Fontname','Arial','Interpreter','latex');
 grid on 
 %MAP LS vs MMSE P_d
 figure(2)
 grid on
 plot(power_noise_ratio_db,(1-(p_md_MAP_LS_final)),'k-s','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,(1-(p_md_MAP_mmse_final)),'r-o','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,(1-(p_md_MAP_ideal_final)),'b-d','LineWidth',2);
 xlabel('$P/N_0$ (dB)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','northwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 grid on 
  %Majority LS vs MMSE P_fa
 figure(3)
 grid on
 plot(power_noise_ratio_db,p_fa_LS_final,'k-s','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,p_fa_mmse_final,'r-o','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,p_fa_ideal_final,'b-d','LineWidth',2);
 xlabel('$P/N_0$ (dB)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','northwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 grid on 
 %Majority LS vs MMSE P_d
 figure(4)
 grid on
 plot(power_noise_ratio_db,(1-(p_md_LS_final)),'k-s','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,(1-(p_md_mmse_final)),'r-o','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,(1-(p_md_ideal_final)),'b-d','LineWidth',2);
 hold on
 xlabel('$P/N_0$ (dB)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','southeast','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 grid on 
 
  %MAP LS vs MMSE P_fa
 figure(5)
 grid on
 plot(power_noise_ratio_db,((p_fa_MAP_LS_final)),'k-s','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,((p_fa_MAP_mmse_final)),'r-o','LineWidth',2);
 hold on
 plot(power_noise_ratio_db,((p_fa_MAP_ideal_final)),'b-d','LineWidth',2);
 xlabel('$P/N_0$ (dB)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','northwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 grid on 
