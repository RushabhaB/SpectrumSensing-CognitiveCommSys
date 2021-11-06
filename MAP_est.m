function [th,CW_p] = MAP_est(th,N0,E_s,h_gain,nSU)
% This function provides the thresholding vector and the probability table
% for MAP combiner
%  Inputs to the function :
%        th : threshold vector
%        N0: Noise power
%        E_s : Energy of the symbols
%        h_gain: Channel power gain
% Outputs the function provides :
%        th : Returns the same th back
%        CW_p :Probability table for the MAP method
gamma = E_s./N0;
%CW_p=zeros(length(th),2,2^(nSU));
for i = 1:length(th)
%Reference: "Cooperative Spectrum Sensing Using Maximum a Posteriori 
%as a Detection Technique for Dynamic Spectrum Access Networks,"

%th(i) = (qfuncinv(fa(i)/2))^2*(N0/2);


for j=1:3
a= sqrt(th(i))/(sqrt(N0(j)/2)*(1+sqrt(gamma(j)))) ;
Pd(j) = gammainc(1,a,'upper');
Pmd(j)=1-Pd(j);
Pfa(j) = 2*qfunc(sqrt(th(i)/N0(j)));
Paf(j) = 1-Pfa(j);
end


CW_p(:,:,i) = [Pd(1)*Pd(2)*Pd(3), Pd(1)*Pd(2)*Pmd(3),Pd(1)*Pmd(2)*Pd(3),Pd(1)*Pmd(2)*Pmd(3),Pmd(1)*Pd(2)*Pd(3),Pmd(1)*Pd(2)*Pmd(3),Pmd(1)*Pmd(2)*Pd(3),Pmd(1)*Pmd(2)*Pmd(3);...
    Pfa(1)*Pfa(2)*Pfa(3), Pfa(1)*Pfa(2)*Paf(3),Pfa(1)*Paf(2)*Pfa(3),Pfa(1)*Paf(2)*Paf(3),Paf(1)*Pfa(2)*Pfa(3),Paf(1)*Pfa(2)*Paf(3),Paf(1)*Paf(2)*Pfa(3),Paf(1)*Paf(2)*Paf(3)]; % Probabilities of each CW. Reference Table 2.


end

end
