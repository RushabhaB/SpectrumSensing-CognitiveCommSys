function [p_md,p_fa] = md_fa_MAP(CW,CW_det,nSamples,nCodeWords)
% This function provides the proability of misdetection and 
% false alarm based on the local FA based on the MAP combiner.
%  Inputs to the function :
%        CW : Actual state of the PU whose size is 1xnCodeWords
%        CW_det : Detected PU state by each SU , matrix size is nSUx nCodeWords
%        nCodeWords : number of codewords
%        nSamples : number of samples for which the channel remains
%        constant
% Outputs the function provides :
%        p_md : Probability of missed detection
%        p_fa : Proability of false alarm

c = CW-CW_det;
md_count = sum(c==1); % Actual is active (1) but estimated is idle (0)
fa_count = sum(c==-1); % Actual is idle (0) but estimated is active (1)

p_md = md_count/(sum(CW==1));
p_fa = fa_count/(sum(CW==0));

end