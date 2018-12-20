function CW = getCW_QLB(qlen)
% Based on (TON 2012, Jiang and Walrand)
rmax = 2;
alpha = 0.23;
%alpha = 0.023;
N_class = 6;
CWmin = 32;
G = log(2);

rarray = rmax - G*(0:1:N_class-1);
[val, id] = min(find(alpha*qlen >= rarray));
CW = CWmin*power(2,val-1);
end