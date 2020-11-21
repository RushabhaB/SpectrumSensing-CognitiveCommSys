function CW_p = MAP_est(th,N0,E_s)

gamma = E_s/N0;
for i = 1:length(th)
a= sqrt(th(i))/(sqrt(N0/2)*(1+sqrt(gamma))) ;

Pmd = gammainc(1,a);
Pd=1-Pmd;

Pfa=2*qfunc(sqrt(th(i)/(N0*0.5)));
Paf = 1-Pfa;

CW_p(i,:,:) = [Pd^3, Pd^2*Pmd,Pd^2*Pmd,Pd^2*Pmd,Pd*Pmd^2,Pd*Pmd^2,Pd*Pmd^2,Pmd^3;Pfa^3,Pfa*Paf^2,...
    Pfa*Paf^2,Pfa*Paf^2,Pfa^2*Paf,Pfa^2*Paf,Pfa^2*Paf,Paf^3];
end
end
