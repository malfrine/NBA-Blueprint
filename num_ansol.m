
tic 

clear 
clc

p0 = sym('p0',[1 12]);
t = sym('t',[1 12]);
v = sym('v', [1 6]);


r = [11.5 11.5 11.5 11 11 11 10.6 10.6 10.5 10.5 10.5 10.5];
p00 =   4.161e+06;
p10 =  -181644;
p01 =  0;%-4.438e+04;
p20 =  37644;
p11 =  0;%      7402;
p02 =  0;%      7950;
K   = 50E6;
C   = 67*30;
a   = 0.01;

%% Developing Objective Equation - Wins

for i = 1:length(p0)
    d(i) = 1/(a*t(i)+1);
    p(i) = p0(i)*d(i);
    pav(i) = int(p(i),0,t(i))/t(i);
end

vars = [p0,t];

for i = 1:length(p0)
    WW(i) = 82*pav(i)/C;
    W(i)  = 82*(pav(i)*t(i)/C - r(i)*t(i)/C);
end

tW = simplify(sum(W));

%% Developing Constraint Equations - Salary and Time

%salary
for i = 1:length(p)
    S(i) = ( p00 + p10*pav(i) + p01*t(i)...
            + p20*pav(i)^2 + p11*pav(i)*t(i) + p02*t(i)^2);
end

tS = simplify(sum(S));

%time
T_PF = t(1) + t(2) + t(3);
T_PG = t(4) + t(5) + t(6);
T_C  = t(7) + t(8);
T_SF = t(9) + t(10);
T_SG = t(11) + t(12);

%% Create Jacobian

jacW     = simplify(jacobian(W,vars));
jacS     = simplify(jacobian(S,vars));
gradT_PF = simplify(jacobian(T_PF,vars));
gradT_PG = simplify(jacobian(T_PG,vars));
gradT_C  = simplify(jacobian(T_C,vars));
gradT_SF = simplify(jacobian(T_SF,vars));
gradT_SG = simplify(jacobian(T_SG,vars));

for i = 1:12
    gradW(1,i) = jacW(i,i);
    gradS(1,i) = jacS(i,i);

end

for i = 13:24
    gradW(1,i) = jacW(i-12,i);
    gradS(1,i) = jacS(i-12,i);
end

%% Develop System of Equations

eqnsObj = gradW - v(1)*gradS - v(2)*gradT_PF - v(3)*gradT_PG ... 
    - v(4)*gradT_C - v(5)*gradT_SF - v(6)*gradT_SG;

eqnsCon1 = tS - K;
eqnsCon2 = T_PF - 48;
eqnsCon3 = T_PG - 48;
eqnsCon4 = T_C - 48;
eqnsCon5 = T_SF - 48;
eqnsCon6 = T_SG - 48;

toc

tic
allEqns = [eqnsObj eqnsCon1 eqnsCon2 eqnsCon3 eqnsCon4 eqnsCon5 eqnsCon6];
vars1 = [vars v];
init = [15*ones(1,12),15*ones(1,12),0.01*ones(1,6)];
sol = vpasolve(allEqns,vars1);
toc


p_sol = [sol.p01 sol.p02 sol.p03 sol.p04 sol.p05 sol.p06 sol.p07 sol.p08...
    sol.p09 sol.p010 sol.p011 sol.p012];
t_sol = [sol.t1 sol.t2 sol.t3 sol.t4 sol.t5 sol.t6 sol.t7 sol.t8 sol.t9 ...
    sol.t10 sol.t11 sol.t12];

%Subsitution

for i = 1:12
    W_sol(i) = subs(W(i),...
    {p0(i) t(i)},...
    {p_sol(i) t_sol(i)});
    S_sol(i) = subs(S(i),...
    {p0(i) t(i)},...
    {p_sol(i) t_sol(i)});
    pav_sol(i) = subs(pav(i),...
        {p0(i) t(i)},...
        {p_sol(i) t_sol(i)});
end

tW_sol = sum(W_sol);
tS_sol = sum(S_sol);
tW_sol






