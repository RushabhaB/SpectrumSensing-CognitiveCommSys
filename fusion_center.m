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

nSU = 3; % Number of Secondary Users (Keep this as 3 to get meaning ful MAP results)
M = 2; %BPSK modulation
%snr = 1:1:25; % Use this array to generate ber vs snr plot

%for i = 1:1:25 %dB

%System Config
N0 = 1;  % Noise power
nCodeWords = 500;
%E_s =100;
%E_s = (10^(i/10))*N0; % Signal energy
%snr = 10*log10(E_s/N0); %SNR of the signal
nSamples = 4 ; % Number of samples across which we assume the H_coeff to be constant
m_p = ones([nSU 1]); % Pilot bits
%CW = zeros([nSU,nCodeWords]); % Actual Codewords received
%CW_det = zeros([nSU,nCodeWords]);% Detected codewords

%power_noise_ratio_db= 0:2:50;
%power_noise_ratio = 10.^(power_noise_ratio_db/10);
E_s = 100;

fa =0.005:0.005:0.05;
th = linspace(10^-4,6,length(fa));
CW_p = MAP_est(th,N0,E_s);
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

for k=1:length(fa)
fa(k)
[x,x_det] = stage1_ED (nSU,nCodeWords,nSamples,E_s,fa(k));
CW = x;
CW_ideal = ones(size(CW));
CW_LS = ones(size(CW));
CW_mmse = ones(size(CW));
CW_detSU = x_det';

pilot_loc = 1:nSamples+1:nCodeWords ;
for i = 1:nSamples+1:nCodeWords
    CW(i) = m_p(1);%The pilots are inserted here as well, though of no consequence it just makes it easier to compare and avoids confusion
    CW_detSU(:,i) = m_p;  %Inserting the pilot vectors into codeword by forcing  the every nSample'th bits to be the pilot
end
%for i 
% Rayleigh Fading Coefficients
H=(randn([nSU nCodeWords])+1j*randn([nSU nCodeWords]))/sqrt(2); % Channel Matrix

% Gaussian noise 
W=sqrt(N0)*(randn([nSU nCodeWords])+randn([nSU nCodeWords])*1j)/sqrt(2);% Noise vector of CSCG noise for SU

for i = 1: nSamples+1 :nCodeWords
x_p = sqrt(E_s).*exp(-1i*pi*2*(m_p)/M); % BPSK Pilot symbol for all the SU

X_p = diag(x_p); % Generating the symbol matrix with the diagonal elements as the data symbols

Y_p = X_p*H(:,i) + W(:,i); % Output symbols at the Fusion Centre (FC)

%Channel  Estimation
H_LS_p(:,ceil(i/(nSamples+1))) = inv(X_p*X_p')*X_p'*Y_p; %MISO Least Square detection using pilots
u = (1./(E_s + N0)).*x_p; %E[h^2] is 2 sigma^2 i.e 2*1/2 = 1. Also ||x||^2 is E_s
H_mmse_p(:,ceil(i/(nSamples+1))) = conj(u).*Y_p; %Standard formula
%dif = H_mmse_p - H_LS_p %Just to check difference... max difference is 0.01 per dimension
end

% Begin Interpolation
for n = 1:nSU
if pilot_loc(1)>1
  slope = (H_LS_p(n,2)-H_LS_p(n,1))/(pilot_loc(2)-pilot_loc(1));
  x = [H_LS_p(n,1)-slope*(pilot_loc(1)-1)  H_LS_p(n,:)]; y = [1 pilot_loc];
  z=  [H_mmse_p(n,1)-slope*(pilot_loc(1)-1)  H_mmse_p(n,:)];
end
if pilot_loc(end)< nCodeWords
  slope = (H_LS_p(n,end)-H_LS_p(n,end-1))/(pilot_loc(end)-pilot_loc(end-1));  
  x = [H_LS_p(n,:)  H_LS_p(n,end)+slope*(nCodeWords-pilot_loc(end))]; y = [pilot_loc nCodeWords];
  z=[H_mmse_p(n,:)  H_mmse_p(n,end)+slope*(nCodeWords-pilot_loc(end))];
end

 H_LS_ipl(n,:) = interp1(y(1:101),x,[1:nCodeWords],'spline');
 H_mmse_ipl(n,:) = interp1(y(1:101),z,[1:nCodeWords],'spline');
end
% End interpolation block

% Begin using the interpolate values to fetermine the data 
for i = 1:nSamples+1:nCodeWords
for j = i+1:i+nSamples

% Data symbols
m = CW_detSU(:,j); % Transmitting the received PU signal by the SU to the FC
xb = sqrt(E_s).*exp(-1i*pi*2*(m)/M); % BPSK Data symbol for all the SU
X_b = diag(xb); % Generating the symbol matrix with the diagonal elements as the data symbols
Y_b = X_b*H(:,j) + W(:,j); % Output symbols at the Fusion Centre (FC)

% Detection using estimated channel 

% Detecting the received bits from the ZF equalisation (LS) of received vector
xb_det = inv(diag(H_LS_ipl(:,j)))*Y_b;
m_det = bpsk_demod(xb_det);

% Detecting the received bits from the ZF equalisation (MMSE) of received vector
xb_det_mmse = inv(diag(H_mmse_ipl(:,j)))*Y_b;
m_det_mmse = bpsk_demod(xb_det_mmse);

%Begin the MAP estimate 

for p = 1:2^nSU
    Ka = sum(abs(Y_b-X_b*H(:,j)).^2) - N0 *log(CW_p(k,1,p));
    Ki = sum(abs(Y_b-X_b*H(:,j)).^2) - N0 *log(CW_p(k,2,p));
    
    Ka_LS = sum(abs(Y_b-X_b*H_LS_ipl(:,j)).^2) - N0 *log(CW_p(k,1,p));
    Ki_LS = sum(abs(Y_b-X_b*H_LS_ipl(:,j)).^2) - N0 *log(CW_p(k,2,p));
    
    Ka_mmse = sum(abs(Y_b-X_b*H_mmse_ipl(:,j)).^2) - N0 *log(CW_p(k,1,p));
    Ki_mmse = sum(abs(Y_b-X_b*H_mmse_ipl(:,j)).^2) - N0 *log(CW_p(k,2,p));
end
Ka_min = min(Ka);
Ki_min = min(Ki);

Ka_min_LS = min(Ka_LS);
Ki_min_LS = min(Ki_LS);

Ka_min_mmse = min(Ka_LS);
Ki_min_mmse = min(Ki_LS);


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
%BER
%[nErr(floor(i/1)),ratio(floor(i/1))] = biterr(m,m_det); %Generating the
%number of errors array and ratio array ( include in snr for loop)
[nErr(j),ratio(j)] = biterr(m,m_det); %Generating the number of errors array and ratio array for ls (without for loop)
[nErr_mmse(j),ratio_mmse(j)] = biterr(m,m_det_mmse); %Generating the number of errors array and ratio array for mmse (without for loop)

CW_detFC(:,j)= m_det;
CW_detFC_mmse(:,j)= m_det_mmse;
end
end

[p_md,p_fa] = md_fa(CW,CW_detFC,nSamples,nCodeWords);
[p_md_mmse,p_fa_mmse] = md_fa(CW,CW_detFC_mmse,nSamples,nCodeWords);

%MAP estimation of Pmd and Pfa

[p_md_MAP_ideal,p_fa_MAP_ideal] = md_fa_MAP(CW,CW_ideal,nSamples,nCodeWords);
[p_md_MAP_LS,p_fa_MAP_LS] = md_fa_MAP(CW,CW_LS,nSamples,nCodeWords);

[p_md_MAP_mmse,p_fa_MAP_mmse] = md_fa_MAP(CW,CW_mmse,nSamples,nCodeWords);

% end MAP estimation

p_md_arr = [p_md_arr p_md];
p_fa_arr = [p_fa_arr p_fa];
p_md_mmse_arr = [p_md_mmse_arr p_md_mmse];
p_fa_mmse_arr = [p_fa_mmse_arr p_fa_mmse];

% MAP figures
p_md_MAP_ideal_arr = [p_md_MAP_ideal_arr p_md_MAP_ideal];
p_fa_MAP_ideal_arr = [p_fa_MAP_ideal_arr p_fa_MAP_ideal];

p_md_MAP_LS_arr = [p_md_MAP_LS_arr p_md_MAP_LS];
p_fa_MAP_LS_arr = [p_fa_MAP_LS_arr p_fa_MAP_LS];

p_md_MAP_mmse_arr = [p_md_MAP_mmse_arr p_md_MAP_mmse];
p_fa_MAP_mmse_arr = [p_fa_MAP_mmse_arr p_fa_MAP_mmse];

end


 figure(1)
 grid on
 hold all
 plot(fa,(p_md_arr),'r-o','LineWidth',2);
 plot(fa,(p_fa_arr),'k-o','LineWidth',2);
 plot(fa,(p_md_mmse_arr),'k--','LineWidth',2);
 plot(fa,(p_fa_mmse_arr),'r--','LineWidth',2);
 xlabel('Local FA probablity','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('Misdetection (LS)', 'False Alarm (LS)','Misdetection (MMSE)', 'False Alarm (MMSE)','Location','southeastoutside','FontSize',12,'Fontname','Arial','Interpreter','latex');

 figure(2)
 grid on
 hold on
 loglog(th,(p_md_MAP_LS_arr),'r-o','LineWidth',2);
 hold on
 loglog(th,(p_md_MAP_mmse_arr),'k--','LineWidth',2);
 hold on
 loglog(th,(p_md_MAP_ideal_arr),'g-.','LineWidth',2);
 
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('Misdetection (LS)','Misdetection (MMSE)','Misdetection (Ideal)','Location','southeastoutside','FontSize',12,'Fontname','Arial','Interpreter','latex');

 figure(3)
 grid on
 hold all
 
 semilogx(th,(p_fa_MAP_LS_arr),'k-o','LineWidth',2);
 
 semilogx(th,(p_fa_MAP_mmse_arr),'r--','LineWidth',2);
 
 semilogx(th,(p_fa_MAP_ideal_arr),'b:','LineWidth',2);
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('False Alarm (LS)','False Alarm (MMSE)','False Alarm (Ideal)','Location','southeastoutside','FontSize',12,'Fontname','Arial','Interpreter','latex');
