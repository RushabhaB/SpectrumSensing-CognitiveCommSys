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

% This file generates the plots for Pfa,Pmd and Pd for MMSE and LS
% estimation methods for each of the combiner techniques i.e. Majority and
% MAP against the local Pfa of the SU which is nothing but the threshold
% value of the ED technique.

% BPSK modulation
bpsk_demod = comm.BPSKDemodulator;% Creating the bpsk demodulator object
nSU = 3;                          % Number of Secondary Users 
M = 2;                            % BPSK modulation

% Local SU System Config
E_s = 100;           % Symbol Energy
N0 = 10;             % Noise power
nCodeWords = 1000;   % Number of CodeWords
nSamples = 4 ;       % Number of samples across which we assume the H_coeff to be constant
m_p = ones([nSU 1]); % Pilot bits
fa = linspace(5.3201e-04,0.9887,25); % Local Pfa vector based on th given in the paper

% Generating the threshold vector based on Pfa and the proability map of
% the MAP scheme based on nSUs 
[th,CW_p] = MAP_est(fa,N0,E_s);

% Number of iterations the entire process runs to average out the results
iter = 5; % 5 for quick testing of the main but we ran it for 5000 iter

for r=1:iter
r  % To display in the command window which iteration it is on

% Initialising all the storage vectors

% Majority combiner arrays
p_md_arr = [];
p_fa_arr = [];
p_md_mmse_arr = [];
p_fa_mmse_arr = [];

% MAP arrays
p_md_MAP_ideal_arr = [];
p_fa_MAP_ideal_arr = [];
p_md_MAP_LS_arr = [];
p_fa_MAP_LS_arr = [];
p_md_MAP_mmse_arr = [];
p_fa_MAP_mmse_arr = [];
p_md_ideal_arr=[];
p_fa_ideal_arr =[];
% End initialisation

for k=1:length(fa)                                         % For each local Pfa the system is run
[x,x_det] = stage1_ED (nSU,nCodeWords,nSamples,E_s,fa(k)); % Energy detection at each SU
CW = x;                                                    % Actual status of the PU at each time instant (1 x nCodeWords)
CW_detSU = x_det';                                         % Status detected by each SU at each time instant (nSU x nCodeWords)

%Initialising the detected PU status for each time instant for each estimation method
CW_ideal = ones(size(CW));
CW_LS = ones(size(CW));
CW_mmse = ones(size(CW));

pilot_loc = 1:nSamples+1:nCodeWords ; %Storing the pilot locations in a vector

for i = 1:nSamples+1:nCodeWords       %Looping over only the pilot positions
    %Inserting Pilots into the SU packets 
    CW(i) = m_p(1);       %The pilots are inserted into the actual status, 
                          %though of no consequence it just makes it easier to compare and avoids confusion
    CW_detSU(:,i) = m_p;  %Inserting the pilot vectors into codeword by forcing the every nSample bits to be the pilot
end

% Channel Matrix of Rayleigh Fading Coefficients 
% Defined only at the pilot samples since we assume that the channel
% remains constant in the time instance between the pilot samples size =(nSU x len(pilot_loc))
H=(randn([nSU ceil(nCodeWords/(nSamples+1))])+1j*randn([nSU ceil(nCodeWords/(nSamples+1))]))/sqrt(2); 

% Gaussian noise : Noise vector of CSCG noise for each SU fo each time
% instant size = (nSU x nCodeWords)
W=sqrt(N0)*(randn([nSU nCodeWords])+randn([nSU nCodeWords])*1j)/sqrt(2);

% Begin estimation at the instants where the pilots are present
for i = 1: nSamples+1 :nCodeWords
x_p = sqrt(E_s).*exp(-1i*pi*2*(m_p)/M); % BPSK Pilot symbol for all the SU

X_p = diag(x_p); % Generating the symbol matrix with the diagonal elements as the data symbols

Y_p = X_p*H(:,ceil(i/(nSamples+1))) + W(:,i); % Output symbols at the Fusion Centre (FC)

%Channel  Estimation
% LS
H_LS_p(:,ceil(i/(nSamples+1))) = inv(X_p*X_p')*X_p'*Y_p; %MISO Least Square detection using pilots

% MMSE
u = (1./(E_s + N0)).*x_p; %E[h^2] is 2 sigma^2 i.e 2*1/2 = 1. Also ||x||^2 is E_s
H_mmse_p(:,ceil(i/(nSamples+1))) = conj(u).*Y_p; %Standard formula
end

% Use the estimated channel at each packet to determine the data symbols
% For all the SUs
for i = 1:nSamples+1:nCodeWords
for j = i+1:i+nSamples

% Data symbols
m = CW_detSU(:,j);                   % Transmitting the data from SU to the FC
xb = sqrt(E_s).*exp(-1i*pi*2*(m)/M); % BPSK Data symbol for all the SU
X_b = diag(xb);                      % Generating the symbol matrix with the diagonal elements as the data symbols
Y_b = X_b*H(:,ceil(i/(nSamples+1))) + W(:,j); % Output CodeWord at the Fusion Centre (FC)

% Using Ideal Channel
% Detecting the received bits from the ZF equalisation  of received vector
xb_det_ideal = inv(diag(H(:,ceil(i/(nSamples+1)))))*Y_b;
% Hard thrsholding based demodulation
m_det_ideal =  bpsk_demod(xb_det_ideal); % Detected status of PU using ideal channel (1xnSU)

% Using LS estimated channel
xb_det = inv(diag(H_LS_p(:,ceil(i/(nSamples+1)))))*Y_b;
m_det = bpsk_demod(xb_det);

% Using MMSE estimated channel
xb_det_mmse = inv(diag(H_mmse_p(:,ceil(i/(nSamples+1)))))*Y_b;
m_det_mmse = bpsk_demod(xb_det_mmse);

%Begin the MAP combiner algorithm
% Ideal Channel combiner
Ka_min = min(sum(abs(Y_b-X_b*H(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,1,:)));
Ki_min= min(sum(abs(Y_b-X_b*H(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,2,:)));

% LS channel combiner
Ka_min_LS = min(sum(abs(Y_b-X_b*H_LS_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,1,:)));
Ki_min_LS = min(sum(abs(Y_b-X_b*H_LS_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,2,:)));

% MMSE channel combiner
Ka_min_mmse = min(sum(abs(Y_b-X_b*H_mmse_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,1,:)));
Ki_min_mmse = min(sum(abs(Y_b-X_b*H_mmse_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,2,:)));

% Checking the PU status
if (Ka_min < Ki_min)
    CW_ideal(j) = 1;
else
    CW_ideal(j) = 0;
end
    
if (Ka_min_LS < Ki_min_LS)
    CW_LS(j) = 1;
else
    CW_LS(j) = 0;
end

if (Ka_min_mmse < Ki_min_mmse)
    CW_mmse(j) = 1;
else
    CW_mmse(j) = 0;
end
%End MAP estimate algorithm

%BER
[nErr(j),ratio(j)] = biterr(m,m_det); %Generating the number of errors array and ratio array for ls (without for loop)
[nErr_mmse(j),ratio_mmse(j)] = biterr(m,m_det_mmse); %Generating the number of errors array and ratio array for mmse (without for loop)

% Detected CodeWord at FC for different channel matrix
CW_detFC_ideal(:,j) = m_det_ideal; % Ideal
CW_detFC(:,j)= m_det;              % LS
CW_detFC_mmse(:,j)= m_det_mmse;    % MMSE
end
end

% Generating the Pfa and Pmd for majority combiner based on the detected CW 
[p_md_ideal,p_fa_ideal] = md_fa(CW,CW_detFC_ideal,nSamples,nCodeWords);
[p_md,p_fa] = md_fa(CW,CW_detFC,nSamples,nCodeWords);
[p_md_mmse,p_fa_mmse] = md_fa(CW,CW_detFC_mmse,nSamples,nCodeWords);

% Pmd and Pfa of MAP combiner  

[p_md_MAP_ideal,p_fa_MAP_ideal] = md_fa_MAP(CW,CW_ideal,nSamples,nCodeWords);
[p_md_MAP_LS,p_fa_MAP_LS] = md_fa_MAP(CW,CW_LS,nSamples,nCodeWords);
[p_md_MAP_mmse,p_fa_MAP_mmse] = md_fa_MAP(CW,CW_mmse,nSamples,nCodeWords);

% End MAP estimation

% Store the Pmd Pfa for each local Pfa for the majority combiner :
% Ideal Channel
p_md_ideal_arr = [p_md_ideal_arr p_md_ideal];
p_fa_ideal_arr = [p_fa_ideal_arr p_fa_ideal];

% LS 
p_md_arr = [p_md_arr p_md];
p_fa_arr = [p_fa_arr p_fa];

% MMSE
p_md_mmse_arr = [p_md_mmse_arr p_md_mmse];
p_fa_mmse_arr = [p_fa_mmse_arr p_fa_mmse];


% Store the Pmd Pfa for each local Pfa for the MAP combiner :
% Ideal Channel
p_md_MAP_ideal_arr = [p_md_MAP_ideal_arr p_md_MAP_ideal];
p_fa_MAP_ideal_arr = [p_fa_MAP_ideal_arr p_fa_MAP_ideal];

% LS channel
p_md_MAP_LS_arr = [p_md_MAP_LS_arr p_md_MAP_LS];
p_fa_MAP_LS_arr = [p_fa_MAP_LS_arr p_fa_MAP_LS];

% MMSE channel
p_md_MAP_mmse_arr = [p_md_MAP_mmse_arr p_md_MAP_mmse];
p_fa_MAP_mmse_arr = [p_fa_MAP_mmse_arr p_fa_MAP_mmse];

end

% Storing the arrays for each iteration

% Majority Combiner : Different channel estimators
% Ideal
p_md_ideal_array(r,:) = p_md_ideal_arr;
p_fa_ideal_array(r,:) = p_fa_ideal_arr;

% LS
p_md_LS_array(r,:) = p_md_arr;
p_fa_LS_array(r,:) = p_fa_arr;

% MMSE
p_md_mmse_array(r,:) = p_md_mmse_arr;
p_fa_mmse_array(r,:) = p_fa_mmse_arr;

% MAP Combiner : Different channel estimators
% Ideal
p_md_MAP_ideal_array(r,:) = p_md_MAP_ideal_arr;
p_fa_MAP_ideal_array(r,:) = p_fa_MAP_ideal_arr;

% LS
p_md_MAP_LS_array(r,:) = p_md_MAP_LS_arr;
p_fa_MAP_LS_array(r,:) = p_fa_MAP_LS_arr;

% MMSE
p_md_MAP_mmse_array(r,:) = p_md_MAP_mmse_arr;
p_fa_MAP_mmse_array(r,:) = p_fa_MAP_mmse_arr;

end

% Taking the mean value across all the iterations
p_md_ideal_final = mean(p_md_ideal_array);
p_fa_ideal_final = mean(p_fa_ideal_array);

p_md_LS_final = mean(p_md_LS_array);
p_fa_LS_final = mean(p_fa_LS_array);

p_md_mmse_final = mean(p_md_LS_array);
p_fa_mmse_final = mean(p_fa_LS_array);

p_md_MAP_ideal_final = mean(p_md_MAP_ideal_array);
p_fa_MAP_ideal_final = mean(p_fa_MAP_ideal_array);

p_md_MAP_LS_final = mean(p_md_MAP_LS_array);
p_fa_MAP_LS_final = mean(p_fa_MAP_LS_array);

p_md_MAP_mmse_final = mean(p_md_MAP_mmse_array);
p_fa_MAP_mmse_final = mean(p_fa_MAP_mmse_array);

% PLOTTING FIGURES

% MAP vs Majority Ideal
 figure(1)
 grid on 
 semilogx(flip(th),flip(1-(p_md_ideal_final)),'k--','LineWidth',2);
 hold on
 semilogx(flip(th),flip(1-(p_md_MAP_ideal_final)),'k-','LineWidth',2);
 hold on
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('Majority Combiner', 'MAP Combiner','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');

 %MAP LS vs MMSE P_d
 figure(2)
 grid on
 semilogx(flip(th),flip(1-(p_md_MAP_LS_final)),'k-s','LineWidth',2);
 hold on
 semilogx(flip(th),flip(1-(p_md_MAP_mmse_final)),'r-d','LineWidth',2);
 hold on
 semilogx(flip(th),flip(1-(p_md_MAP_ideal_final)),'b-+','LineWidth',2);
 hold on
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')

 %MAP LS vs MMSE P_fa
 figure(3)
 grid on
 semilogx(flip(th),flip((p_fa_MAP_LS_final)),'k-s','LineWidth',2);
 hold on
 semilogx(flip(th),flip((p_fa_MAP_mmse_final)),'r-+','LineWidth',2);
 hold on
 semilogx(flip(th),flip((p_fa_MAP_ideal_final)),'b-d','LineWidth',2);
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 
 
 %Majority LS vs MMSE P_d
 figure(4)
 grid on
 semilogx(flip(th),flip(1-(p_md_LS_final)),'k-s','LineWidth',2);
 hold on
 semilogx(flip(th),flip(1-(p_md_mmse_final)),'r-+','LineWidth',2);
 hold on
 semilogx(flip(th),flip(1-(p_md_ideal_final)),'b-d','LineWidth',2);
 hold on
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')

 %Majority LS vs MMSE P_fa
 figure(5)
 grid on
 semilogx(flip(th),flip((p_fa_LS_final)),'k-s','LineWidth',2);
 hold on
 semilogx(flip(th),flip((p_fa_mmse_final)),'r-+','LineWidth',2);
 hold on 
 semilogx(flip(th),flip((p_fa_ideal_final)),'b-d','LineWidth',2);
 hold on
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 

 
 % Majority MD vs FA
 figure(6)
 grid on 
 hold all
 plot(fa,p_md_LS_final,'k-d','LineWidth',2);
 hold on
 plot(fa,p_md_mmse_final,'r-+','LineWidth',2);
 hold on
 plot(fa,p_md_ideal_final,'b-d','LineWidth',2);
 hold on
 xlabel('Probablity of local SU false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of misdetection ($P_{md}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','northwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 

 % MAP MD vs FA
 figure(7)
 grid on 
 hold all
 plot(fa,p_md_MAP_LS_final,'k-o','LineWidth',2);
 hold on
 plot(fa,p_md_MAP_mmse_final,'r-+','LineWidth',2);
 hold on
 plot(fa,p_md_MAP_ideal_final,'b-d','LineWidth',2);
 hold on
 xlabel('Probablity of local SU false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of misdetection ($P_{md}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Ideal','Location','northeast','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 

 
 
%  figure(5)
%  grid on
%  sgtitle('False Alarm Probability : MAP combiner')
%  subplot(311)
%  semilogx(th,(p_fa_MAP_LS_final),'k-o','LineWidth',2);
%  legend('False Alarm (LS)','Location','best','FontSize',10,'Fontname','Arial','Interpreter','latex')
%  subplot(312)
%  semilogx(th,(p_fa_MAP_mmse_final),'r--','LineWidth',2);
%  legend('False Alarm (MMSE)','Location','best','FontSize',10,'Fontname','Arial','Interpreter','latex')
%  subplot(313)
%  semilogx(th,(p_fa_MAP_ideal_final),'b:','LineWidth',2);
%  legend('False Alarm (Ideal)','Location','best','FontSize',10,'Fontname','Arial','Interpreter','latex')
%  xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
%  ylabel('Probablity','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')

