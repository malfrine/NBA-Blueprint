clc
clear

p = sym('p',[1 12]);
t = sym('t',[1 12]);
r = sym('r',[1 12]);
v = sym('v', [1 6]);
vars = [p,t];

syms as bs cs K a C


%% Developing Objective Equation - Wins

for i = 1:length(p) 
    W(i) = p(i)*log(1+a*t(i))/a/C - r(i)*t(i)/C;
end

tW = simplify(sum(W));

%% Developing Contraint Equations - Salary and Time

%salary
for i = 1:length(p)
    S(i) = as*p(i).^2+bs*p(i)+cs;  
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

eqnsCon1 = S - K;
eqnsCon2 = T_PF - 48;
eqnsCon3 = T_PG - 48;
eqnsCon4 = T_C - 48;
eqnsCon5 = T_SF - 48;
eqnsCon6 = T_SG - 48;

allEqns = [eqnsObj eqnsCon1 eqnsCon2 eqnsCon3 eqnsCon3 eqnsCon4 ...
    eqnsCon5 eqnsCon6];

% sol = solve(allEqns,vars);

%% Sensitivity Analysis

dtWda = diff(tW,a);

Vars = {C a p(1) p(2) p(3) p(4) p(5) p(6) p(7) p(8) p(9) p(10) p(11) p(12) ...
    t(1) t(2) t(3) t(4) t(5) t(6) t(7) t(8) t(9) t(10) t(11) t(12)};

Sub_Vars{1} = {2E7 0.01 12.897	12.839	12.897	12.896	12.898	12.838	12.410	12.411	12.410	12.411	12.411	12.410...
16.032	15.936	16.032	16.032	16.032	15.936	24.000	24.000	24.000	24.000	24.000	24.000};

Sub_Vars{2} = {2E7 0.01 15.372	15.372	15.295	15.372	15.372	15.295	14.791	14.791	14.791	14.791	14.791	14.791...
16.033	16.033	15.933	16.033	16.033	15.933	24.000	24.000	24.000	24.000	24.000	24.000};

Sub_Vars{3} = {2E7 0.01 17.421	17.339	17.421	17.421	17.421	17.339	16.762	16.762	16.762	16.762	16.762	16.762...
16.034	15.933	16.034	16.034	16.034	15.933	24.000	24.000	24.000	24.000	24.000	24.000};

W_0 = [17.40219683 41.05723232 60.65159703];

for i = 1:3
ratio = double(subs(dtWda,Vars,Sub_Vars{i}));
Sp(i) = double(W_0(i)/0.01*ratio);
end



save 'NBA_Function'











