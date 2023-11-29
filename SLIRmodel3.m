%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SLIR function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIRmodel combines all the SLIR equations into one system.  
% the dependent variables are stored in y as:
% y(1) = S (susceptible population)
% y(2) = L (latent population)
% y(3) = I (infectious population)
% y(4) = R (recovered/removed population)
%
% with parameters
% p(1) = beta (rate of new infections)
% p(2) = mu_l (inverse length of latent period in days)
% p(3) = mu_i (inverse length of the infectious period in days)
% p(4) = k (population growth rate)
% p(5) = e (rate of new infections from external sources)

function [dydt] = SLIRmodel3(t,y,p,beta)
    %assign parameter
    mu_l = p(1,1);
    mu_i = p(1,2);
    e    = p(1,3);
    load EnvironmentalForcing.mat;
    Ap = p(1,4);
    day = p(1,5);
    TE = -0.35968 + 0.10789*T(t) - 0.00214.*T(t)^2;
    %assign variables
    S = y(1);
    L = y(2);
    I = y(3);
    R = y(4);
    Pb = y(5);

    dydt(6) = (0.1724*Pb-.0000212*Pb^2)*TE; %dPb/dt
    dydt(5) = (1.33*day)*TE;%dPl/dt
    dydt(1) = -beta*S*I+(dydt(6)+dydt(5))*(1/Ap); %dS/dt
    dydt(2) = S*I-mu_l*L+e;%dL/dt
    dydt(3) = mu_l*L-mu_i^-1*I;%dI/dt
    dydt(4) = mu_i^-1*I;%dR/dt

end