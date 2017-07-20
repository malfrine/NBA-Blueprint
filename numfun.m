function [c,ceq] = confuneq(var)
p00 =  4E6;
p10 =  -181644;%-5.226e+05;%
p01 =  0;%-4.438e+04;
p20 =  30167;%2.248e+04;%
p11 =  0;%      7402;%
p02 =  0;%      7950;%
K   =  130E6;
a   =  0.01;
% Nonlinear equality constraints
for i = 1:12
    poi    = var(i);
    ti     = var(i+12);
    pavi   = poi*(1-exp(-a*ti))/a/ti;
    ceqi(i) = ( p00 + p10*pavi + p01*ti...
            + p20*pavi^2 + p11*pavi*ti + p02*ti^2);
    c(i) = var(i);
    c(i+12) = var(i+12);
end

ceq(1) = var(1+12) + var(2+12) + var(3+12) - 48;
ceq(2) = var(4+12) + var(5+12) + var(6+12) - 48;
ceq(3) = var(7+12) + var(8+12) - 48;
ceq(4) = var(9+12) + var(10+12) - 48;
ceq(5) = var(11+12) + var(12+12) - 48;
ceq(6) = sum(ceqi) - K;


end