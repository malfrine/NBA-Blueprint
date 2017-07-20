clear
clc

tic

var0 = [15*ones(1,12) 16*ones(1,12)]; % Make a starting guess at the solution
   opts = optimset('Display','iter','Algorithm','interior-point',...
                  'MaxFunEval',inf,'MaxIter',100);
[x,fval] = fmincon(@objfun,var0,[],[],[],[],[],[],... 
   @confuneq,opts);

toc

export = [x(1:12);x(13:24)];