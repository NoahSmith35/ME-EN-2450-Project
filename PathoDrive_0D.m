% driver script for our zero dimensional pathogen simulation
%
% set simulation constants
beta = 2;    %rate infection increases (1/day)
mu_L = 0.5;  %rate latent period ends (inverse of number of days latent)
mu_I = 0.25; %rate infection clears (inverse of number of days infectious)
days = 60;   %length of simulation (day)
dt   = 0.05; %timestep (fraction of a day)
S_i  = 1;    %initial size of the population (normalized)
L_i  = S_i*0.001; %initial fraction of population that is latent
I_i  = 0;    %initial fraction of population that is infectious
R_i  = I_i*mu_I;  %initial fraction of population that is recovered
k    = 0.01; %population growth rate (fraction per day)
e    = 0.005;%rate of introduction from external sources

% call the pathogen function
tic
[S,L,I,R,time]=PathogenGrowth_0D(S_i,L_i,I_i,R_i,beta,mu_L,mu_I,k,e,days,dt);
toc

% plot results
FSize = 14; %fontsize for plots
figure;plot(time,S+L+I+R,'-k',time,S,'-.m',time,L,'--g',time,I,':b',...
    time,R,'-.r','LineWidth',2);
legend({'Total Population';'Susceptible';'Latent';'Infected';'Removed'},'Location','NorthWest');
xlabel('time (days)','Fontsize',FSize);
ylabel('Population (fraction of initial)','Fontsize',FSize)
set(gca,'Fontsize',FSize);
box on;grid on
