function tW = objfun(var)
r = [11.5 11.5 11.5 11 11 11 10.6 10.6 10.5 10.5 10.5 10.5];
C = 67*30;
a = 0.01;

for i = 1:12
    poi  = var(i);
    ti   = var(i+12);
    d(i) = exp(-a*ti);
    pti  = poi*d(i);
    pav  = poi*(1-exp(-a*ti))/a/ti;
    W(i) = 82*(pav*ti - r(i)*ti) / C;
end
tW = -sum(W) ;
end