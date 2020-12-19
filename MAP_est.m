function [th,CW_p] = MAP_est(fa,N0,E_s,h_gain,nSU)
% This function provides the thresholding vector and the probability table
% for MAP combiner
%  Inputs to the function :
%        fa : False alarm vector
%        N0: Noise power
%        E_s : Energy of the symbols
%        h_gain: Channel power gain
% Outputs the function provides :
%        th : Threshold for ED method based on the FA
%        CW_p :Probability table for the MAP method
gamma = E_s/N0;
CW_p=ones(length(fa),2,2^(nSU));
for i = 1:length(fa)
%Reference: "Cooperative Spectrum Sensing Using Maximum a Posteriori 
%as a Detection Technique for Dynamic Spectrum Access Networks,"

th(i) = (qfuncinv(fa(i)/2))^2*(N0/2);
a= sqrt(th(i))/(sqrt(N0/2)*(1+sqrt(gamma))) ;

Pd = gammainc(1,a,'upper');
Pmd=1-Pd;

Pfa=fa(i);
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
