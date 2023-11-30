clear
clc
load EnvironmentalForcing.mat;
beta_max = 1;
mu_L_min = 6;
mu_I = 10;
e = 0.001;
Ap = 5000;
P_i = 1.33 * 30 * (-0.35968 + 0.10789 * 15 - 0.00214 * 15 * 15) * 30;
S_i = P_i./Ap;
L_i = 0.01*S_i;
I_i = 0;
R_i = mu_I*I_i;
day = 61;
B_i = 1;
y = [S_i L_i I_i R_i, 1];


T_beta = zeros(1,length(T));
for i = 1:length(T)
  t = T(i);
  if t<35
T_beta(i) = 0.000241.*t.^2.06737 .* (35-t).^0.72859;
% extravar = T_beta/24;
   end
end
j = 1;
mu_L = zeros(1,length(T));
for i= 1:length(T)
    mu_L(i) = sum(T_beta(j:i));
    while mu_L(i)>mu_L_min
        j = j+1;
        mu_L(i)= sum(T_beta(j:i));
    end
end
beta = beta_max.* T_beta;

for i = 1:length(T)
    p = [mu_L(i),mu_I,e,Ap,day];
    for j = 1:4
    ansArray(j,i)=rk4(@dydt,tspan, y(j));
    end
end
