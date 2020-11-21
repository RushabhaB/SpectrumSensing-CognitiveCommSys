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
N0 = 10;  % Noise power
nCodeWords = 1000;
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

fa = linspace(5.3201e-04,0.9887,25); % FA vector values taken from the paper based on th : 10^-4 - 6 W
%th = linspace(10^-4,6,length(fa));
[th,CW_p] = MAP_est(fa,N0,E_s);


iter = 5;

for r=1:iter
    r
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
for k=1:length(fa)
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
H=(randn([nSU ceil(nCodeWords/(nSamples+1))])+1j*randn([nSU ceil(nCodeWords/(nSamples+1))]))/sqrt(2); % Channel Matrix

% Gaussian noise 
W=sqrt(N0)*(randn([nSU nCodeWords])+randn([nSU nCodeWords])*1j)/sqrt(2);% Noise vector of CSCG noise for SU

for i = 1: nSamples+1 :nCodeWords
x_p = sqrt(E_s).*exp(-1i*pi*2*(m_p)/M); % BPSK Pilot symbol for all the SU

X_p = diag(x_p); % Generating the symbol matrix with the diagonal elements as the data symbols

Y_p = X_p*H(:,ceil(i/(nSamples+1))) + W(:,i); % Output symbols at the Fusion Centre (FC)

%Channel  Estimation
H_LS_p(:,ceil(i/(nSamples+1))) = inv(X_p*X_p')*X_p'*Y_p; %MISO Least Square detection using pilots
u = (1./(E_s + N0)).*x_p; %E[h^2] is 2 sigma^2 i.e 2*1/2 = 1. Also ||x||^2 is E_s
H_mmse_p(:,ceil(i/(nSamples+1))) = conj(u).*Y_p; %Standard formula
%dif = H_mmse_p - H_LS_p %Just to check difference... max difference is 0.01 per dimension
end

% % Begin Interpolation
% for n = 1:nSU
% if pilot_loc(1)>1
%   slope = (H_LS_p(n,2)-H_LS_p(n,1))/(pilot_loc(2)-pilot_loc(1));
%   x = [H_LS_p(n,1)-slope*(pilot_loc(1)-1)  H_LS_p(n,:)]; y = [1 pilot_loc];
%   z=  [H_mmse_p(n,1)-slope*(pilot_loc(1)-1)  H_mmse_p(n,:)];
% end
% if pilot_loc(end)< nCodeWords
%   slope = (H_LS_p(n,end)-H_LS_p(n,end-1))/(pilot_loc(end)-pilot_loc(end-1));  
%   x = [H_LS_p(n,:)  H_LS_p(n,end)+slope*(nCodeWords-pilot_loc(end))]; y = [pilot_loc nCodeWords];
%   z=[H_mmse_p(n,:)  H_mmse_p(n,end)+slope*(nCodeWords-pilot_loc(end))];
% end
% 
%  H_LS_ipl(n,:) = interp1(y(1:101),x,[1:nCodeWords],'spline');
%  H_mmse_ipl(n,:) = interp1(y(1:101),z,[1:nCodeWords],'spline');
% end
% % End interpolation block

% Begin using the interpolate values to fetermine the data 
for i = 1:nSamples+1:nCodeWords
for j = i+1:i+nSamples

% Data symbols
m = CW_detSU(:,j); % Transmitting the received PU signal by the SU to the FC
xb = sqrt(E_s).*exp(-1i*pi*2*(m)/M); % BPSK Data symbol for all the SU
X_b = diag(xb); % Generating the symbol matrix with the diagonal elements as the data symbols
Y_b = X_b*H(:,ceil(i/(nSamples+1))) + W(:,j); % Output symbols at the Fusion Centre (FC)

% Detection using estimated channel 

% Detecting the received bits from the ZF equalisation (LS) of received vector
xb_det_ideal = inv(diag(H(:,ceil(i/(nSamples+1)))))*Y_b;
m_det_ideal =  bpsk_demod(xb_det_ideal);

xb_det = inv(diag(H_LS_p(:,ceil(i/(nSamples+1)))))*Y_b;
m_det = bpsk_demod(xb_det);

% Detecting the received bits from the ZF equalisation (MMSE) of received vector
xb_det_mmse = inv(diag(H_mmse_p(:,ceil(i/(nSamples+1)))))*Y_b;
m_det_mmse = bpsk_demod(xb_det_mmse);

%Begin the MAP estimate 
Ka_min = min(sum(abs(Y_b-X_b*H(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,1,:)));
Ki_min= min(sum(abs(Y_b-X_b*H(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,2,:)));

Ka_min_LS = min(sum(abs(Y_b-X_b*H_LS_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,1,:)));
Ki_min_LS = min(sum(abs(Y_b-X_b*H_LS_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,2,:)));

Ka_min_mmse = min(sum(abs(Y_b-X_b*H_mmse_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,1,:)));
Ki_min_mmse = min(sum(abs(Y_b-X_b*H_mmse_p(:,ceil(i/(nSamples+1)))).^2) - N0 *log(CW_p(k,2,:)));

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
%End MAP estimate

%BER

[nErr(j),ratio(j)] = biterr(m,m_det); %Generating the number of errors array and ratio array for ls (without for loop)
[nErr_mmse(j),ratio_mmse(j)] = biterr(m,m_det_mmse); %Generating the number of errors array and ratio array for mmse (without for loop)

CW_detFC_ideal(:,j) = m_det_ideal;
CW_detFC(:,j)= m_det;
CW_detFC_mmse(:,j)= m_det_mmse;
end
end

[p_md_ideal,p_fa_ideal] = md_fa(CW,CW_detFC_ideal,nSamples,nCodeWords);
[p_md,p_fa] = md_fa(CW,CW_detFC,nSamples,nCodeWords);
[p_md_mmse,p_fa_mmse] = md_fa(CW,CW_detFC_mmse,nSamples,nCodeWords);

%MAP estimation of Pmd and Pfa

[p_md_MAP_ideal,p_fa_MAP_ideal] = md_fa_MAP(CW,CW_ideal,nSamples,nCodeWords);
[p_md_MAP_LS,p_fa_MAP_LS] = md_fa_MAP(CW,CW_LS,nSamples,nCodeWords);

[p_md_MAP_mmse,p_fa_MAP_mmse] = md_fa_MAP(CW,CW_mmse,nSamples,nCodeWords);

% end MAP estimation
p_md_ideal_arr = [p_md_ideal_arr p_md_ideal];
p_fa_ideal_arr = [p_fa_ideal_arr p_fa_ideal];

p_md_arr = [p_md_arr p_md];
p_fa_arr = [p_fa_arr p_fa];

p_md_mmse_arr = [p_md_mmse_arr p_md_mmse];
p_fa_mmse_arr = [p_fa_mmse_arr p_fa_mmse];


% % MAP figures
p_md_MAP_ideal_arr = [p_md_MAP_ideal_arr p_md_MAP_ideal];
p_fa_MAP_ideal_arr = [p_fa_MAP_ideal_arr p_fa_MAP_ideal];
% 
p_md_MAP_LS_arr = [p_md_MAP_LS_arr p_md_MAP_LS];
p_fa_MAP_LS_arr = [p_fa_MAP_LS_arr p_fa_MAP_LS];
% 
p_md_MAP_mmse_arr = [p_md_MAP_mmse_arr p_md_MAP_mmse];
p_fa_MAP_mmse_arr = [p_fa_MAP_mmse_arr p_fa_MAP_mmse];



end

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
 %subplot(3,1,1)
 semilogx(flip(th),flip(1-(p_md_MAP_LS_final)),'k-s','LineWidth',2);
 %legend('Detection (LS)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(312)
 hold on
 semilogx(flip(th),flip(1-(p_md_MAP_mmse_final)),'r--','LineWidth',2);
 %legend('Detection (MMSE)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(313)
 %semilogx(flip(th),flip(1-(p_md_MAP_ideal_final)),'b:','LineWidth',2);
 %legend('Detection (Ideal)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %han=axes(fig,'visible','off'); 
 %han.XLabel.Visible='on';
 %han.YLabel.Visible='on';
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')

 %MAP LS vs MMSE P_fa
 figure(3)
 grid on
 %subplot(3,1,1)
 semilogx(flip(th),flip((p_fa_MAP_LS_final)),'k-s','LineWidth',2);
 %legend('Detection (LS)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(312)
 hold on
 semilogx(flip(th),flip((p_fa_MAP_mmse_final)),'r--','LineWidth',2);
 %legend('Detection (MMSE)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(313)
 %semilogx(flip(th),flip(1-(p_md_MAP_ideal_final)),'b:','LineWidth',2);
 %legend('Detection (Ideal)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %han=axes(fig,'visible','off'); 
 %han.XLabel.Visible='on';
 %han.YLabel.Visible='on';
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('MAP combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 
 
 %Majority LS vs MMSE P_d
 figure(4)
 grid on
 %subplot(3,1,1)
 semilogx(flip(th),flip(1-(p_md_LS_final)),'k-s','LineWidth',2);
 %legend('Detection (LS)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(312)
 hold on
 semilogx(flip(th),flip(1-(p_md_mmse_final)),'r--','LineWidth',2);
 %legend('Detection (MMSE)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(313)
 %semilogx(flip(th),flip(1-(p_md_MAP_ideal_final)),'b:','LineWidth',2);
 %legend('Detection (Ideal)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %han=axes(fig,'visible','off'); 
 %han.XLabel.Visible='on';
 %han.YLabel.Visible='on';
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of detection ($P_d$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')

 %Majority LS vs MMSE P_fa
 figure(5)
 grid on
 %subplot(3,1,1)
 semilogx(flip(th),flip((p_fa_LS_final)),'k-s','LineWidth',2);
 %legend('Detection (LS)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(312)
 hold on
 semilogx(flip(th),flip((p_fa_mmse_final)),'r--','LineWidth',2);
 %legend('Detection (MMSE)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %subplot(313)
 %semilogx(flip(th),flip(1-(p_md_MAP_ideal_final)),'b:','LineWidth',2);
 %legend('Detection (Ideal)','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 %han=axes(fig,'visible','off'); 
 %han.XLabel.Visible='on';
 %han.YLabel.Visible='on';
 xlabel('Threshold(W)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Location','southwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 

 
 % Majority MD vs FA
 figure(6)
 grid on 
 hold all
 plot(fa,p_md_LS_final,'k-d','LineWidth',2);
 hold on
 plot(fa,p_md_mmse_final,'r--','LineWidth',2);
 xlabel('Probablity of local SU false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of misdetection ($P_{md}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Location','northwest','FontSize',10,'Fontname','Arial','Interpreter','latex');
 title('Majority combiner','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 

 % MAP MD vs FA
 figure(7)
 grid on 
 hold all
 plot(fa,p_md_MAP_LS_final,'k-d','LineWidth',2);
 hold on
 plot(fa,p_md_MAP_mmse_final,'r--','LineWidth',2);
 xlabel('Probablity of local SU false alarm ($P_{fa}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 ylabel('Probablity of misdetection ($P_{md}$)','FontSize',12,'FontWeight','bold','Color','k','Fontname', 'Arial','Interpreter', 'latex')
 legend('LS Estimate', 'MMSE Estimate','Location','northeast','FontSize',10,'Fontname','Arial','Interpreter','latex');
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

