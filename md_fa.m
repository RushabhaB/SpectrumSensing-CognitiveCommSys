% Function to give P_md (misdetection) and P_fa (false alarm) 

function [p_md,p_fa] = md_fa(CW,CW_det,nSamples,nCodeWords)

%a = sum(CW);
b = sum(CW_det);

%a(a<2) = 0; % Actual PU status - Idle
%a(a>0) = 1; % Actual PU status - Active
a=CW;
b(b<2) = 0; % Estimated PU status - Idle
b(b>0) = 1; % Estimated PU status - Active

pilots_len = length(1: nSamples+1 :nCodeWords); % Number of pilot codewords
actual_data_len = nCodeWords - pilots_len; % Number of actual data codewords

c = a-b;
md_count = sum(c==1); % Actual is active (1) but estimated is idle (0)
fa_count = sum(c==-1); % Actual is idle (0) but estimated is active (1)

p_md = md_count/(actual_data_len);
p_fa = fa_count/(actual_data_len);

end

