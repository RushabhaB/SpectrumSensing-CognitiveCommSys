function [th,CW_p] = MAP_est_Es_change(fa,N0,E_s)
%fa = linspace(5.3201e-04,0.9887,25);
%N0=1;
%E_s = 100;


for i = 1:length(E_s)
gamma = E_s(i)/N0;
th = (qfuncinv(fa/2))^2*(N0/2);
a= sqrt(th)/(sqrt(N0/2)*(1+sqrt(gamma))) ;

Pd = gammainc(1,a,'upper');
Pmd=1-Pd;

Pfa=fa;
Paf = 1-Pfa;

CW_p(i,:,:) = [Pd^3, Pd^2*Pmd,Pd^2*Pmd,Pd^2*Pmd,Pd*Pmd^2,Pd*Pmd^2,Pd*Pmd^2,Pmd^3;Pfa^3,Pfa*Paf^2,...
    Pfa*Paf^2,Pfa*Paf^2,Pfa^2*Paf,Pfa^2*Paf,Pfa^2*Paf,Paf^3];
end
end
