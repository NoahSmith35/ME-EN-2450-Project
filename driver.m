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
day = 
B_i = 1;
T_beta = zeros(1,length(T));
for i = 1:length(T)
  t = T(i);
  if t<35
T_beta(i) = 0.000241.*t.^2.06737 .* (35-t).^0.72859;
  end
end
beta = beta_max * T_beta;
sum1 =sum( T_beta);






