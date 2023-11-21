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

function [dydt] = SLIRmodel(t,y,p)
    %assign parameters
    beta = p(1);
    mu_l = p(2);
    mu_i = p(3);
    k    = p(4);
    e    = p(5);
    %assign variables
    S = y(1);
    L = y(2);
    I = y(3);
    R = y(4);
    
    dydt(1) = -beta*S*I+k;
    dydt(2) = beta*S*I-mu_l*L+e;
    dydt(3) = mu_l*L-mu_i*I;
    dydt(4) = mu_i*I;
end