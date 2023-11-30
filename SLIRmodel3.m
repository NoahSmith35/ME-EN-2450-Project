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

function [outfun] = SLIRmodel3(t,yvect,j,i)
    %assign parameter
    load SLIRPvars.mat
    TE = -0.35968 + 0.10789*T(i) - 0.00214.*T(i)^2;
    
    %assign variables
    S = yvect(1,1);
    L = yvect(2,1);
    I = yvect(3,1);
    R = yvect(4,1);
    Pb = yvect(5,1);

    dPbdt = (0.1724*Pb-.0000212*Pb^2)*TE; %dPb/dt
    dPldt = (1.33*(t/24))*TE;%dPl/dt
    dydt(5) = dPbdt+dPldt;
    dydt(1) = -Beta(i)*S*I+(dydt(5))*(1/Ap); %dS/dt
    dydt(2) = Beta(i)*S*I-(mu_L(i)^-1)*L+e;%dL/dt
    dydt(3) = (mu_L(i)^-1)*L-mu_I^-1*I;%dI/dt
    dydt(4) = mu_I^-1*I;%dR/dt
    outfun = dydt(j);
end
