% Cooperative Cognitive Radio Communication
% MATLAB code for assignment in AET G641 @ BITS Pilani
% Instructor: B. Sainath

% Students:
% Vandana Prasad - 2019H1240092P
% Rushabha Balaji - 2017A3PS0220P
% Vinay U Pai - 2017A3PS0131P

% Function for MAP combining
% Refer MAP_est.m for funtionality , the only change is E_s is being varied
% instead of local Pfa.

function [th,CW_p] = MAP_est_Es_change(fa,N0,E_s)

CW_p=ones(length(fa),2,2^(nSU));

for i = 1:length(E_s)
%Reference: "Cooperative Spectrum Sensing Using Maximum a Posteriori 
%as a Detection Technique for Dynamic Spectrum Access Networks,"

gamma = E_s(i)/N0; %SNR 
th = (qfuncinv(fa/2))^2*(N0/2); % Threshold calculated from P_fa
a= sqrt(th)/(sqrt(N0/2)*(1+sqrt(gamma))) ;
Pd = gammainc(1,a,'upper'); % P_d calculated from Threshold
Pmd=1-Pd; %P_md

Pfa=fa;
Paf = 1-Pfa;

%CW_p(i,:,:) = [Pd^3, Pd^2*Pmd,Pd^2*Pmd,Pd^2*Pmd,Pd*Pmd^2,Pd*Pmd^2,Pd*Pmd^2,Pmd^3;Pfa^3,Pfa*Paf^2,...
%    Pfa*Paf^2,Pfa*Paf^2,Pfa^2*Paf,Pfa^2*Paf,Pfa^2*Paf,Paf^3]; % Probabilities of each CW. Reference Table 2.

CW_p(i,1,1) = Pd^nSU;
CW_p(i,2,1) = Pfa^nSU;
CW_p(i,1,end) = Pmd^nSU;
CW_p(i,2,end) = Paf^nSU;
init = 2;
len1 = 0;

for j=1:(nSU-1)
    len2 = nchoosek(nSU,j)-1;
    init = init +len1;
    CW_p(i,1,init:init+len2)=repmat(Pmd^j*Pd^(nSU-j),nchoosek(nSU,j),1)';
    CW_p(i,2,init:init+len2)=repmat(Paf^j*Pfa^(nSU-j),nchoosek(nSU,j),1)';
    len1 = len2+1;
end

end
end
