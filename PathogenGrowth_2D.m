% function to compute the increase in icidence of a plant pathogen (or any
% pathogen with stationary hosts).  This code is 2-D because it predicts 
% the growth of the pathogen for a 2-D array of hosts each with a 
% population that grows in size with time.  Additionally, it calculates 
% airborne transmission of disease between hosts in the 2D array.  It is 
% a modified version of the PathogenGrowth_0D function used in lab 09.
%
% The primary equations governing this are compound interest equations for
% growth:
%
% dS/dt = -beta*S*I+dP/dt
%
% dL/dt = beta*S*I-mu_L*L+dE/dt
%
% dI/dt = mu_L*L*I-mu_I*I
%
% dP/dt = dB/dt + d(leaf)/dt 
%
% dR/dt = mu_I*I
%
% dE/dt = e
%
% dF/dt = Gamma*exp(alpha*I) - F*R_frac
%
% with:
% S    = susceptible fraction of population (susceptible tissue)
% beta = rate of infection growth for healthy population (fraction per day)
% L    = fraction of tissue infected and latent (e.g., dormant before infection)
% I    = fraction of tissue infected and producing inoculum
% R    = fraction of tissue recovered (or removed) from population
% P    = size of the total population (plant surface area)
% B    = surface area of berries
% E    = amount (fraction of population) of introduced disease from external sources
% F    = size of the spreading population, e.g. sporulating for a fungus 
% t    = time
% mu_L = rate of decrease in latent population (fraction per day)
% mu_I = rate of decrease in infectious population (fraction per day)
% e    = rate of import from external sources
%
% inputs: vine (structure containing the initial size of susceptible population); beta; mu;
% tspan (days to simulate array);
% output: S,L,I,R,P,E,time (vector of simulation times), and B

function [vine,infects,infectsFound,tFound,cost] = PathogenGrowth_2D(vine,beta_max,mu_L_target,mu_I,A,...
    eta,kappa,xi,Gamma,alpha,T,U,V,tspan,onlyScout,scoutingMethod)

%declare global variables
global NpX NpY Nsteps

% set parameters in a cell array
p{1} = beta_max; %(max rate of new infections)
p{2} = 1/mu_I;   %(inverse length of the infectious period in days)
p{3} = T;        %(temperature in C)
p{4} = tspan;    %(time in days array)
p{5} = A;        %(total plant surface area at reference time)
p{6} = sqrt(U.^2+V.^2);  %windspeed
p{7} = atand(V./U);      %wind direction
p{8} = eta;      %release fraction scale factor
p{9} = kappa;    %release fraction scale factor   
p{10}= xi;       %release fraction offset
p{11}= Gamma;    %spore production multiple
p{12}= alpha;    %spore production 2nd factor
cost = 0;        %Intialize cost
findSwitch = 0;  %if the disease has been found
tFound = Inf;    %time the disease was found
infectsFound = 0;
infects = zeros(NpX,NpY); % the loction the disease was found 
% declare function handles
odefun = @(t,y,e,g) SLIRPE_model(t,y,e,g,p);
% loop over timesteps (starting at 2)
for t=2:Nsteps

%     if mod(t-1,24)==0
        disp(['day=',num2str(tspan(t),'%.2f'),' infected plants=',int2str(sum([vine.IsInfect]))])
%     end
    dt=tspan(t)-tspan(t-1); %timestep

    %update list of infected vines
    ActiveVines = find([vine.IsInfect]);
    %initialize deposition flux for timestep
    DepFlux_sum=zeros(1,NpX*NpY);

    for idx = ActiveVines
        %calculate positions of other vines relative to each infected plant
        Xplume = [vine.X]-vine(idx).X;
        Yplume = [vine.Y]-vine(idx).Y;
        
        DepFlux = GaussianPlumeDep(Xplume,Yplume,p{6}(t),p{7}(t),vine(idx).S(t-1),vine(idx).F(t-1));
        NNaNInd = ~isnan(DepFlux);
        DepFlux_sum(NNaNInd) = DepFlux_sum(NNaNInd) + DepFlux(NNaNInd);
    end
    
    %now for all vines
    for i=1:NpX
        for j=1:NpY
            cnt=i+(j-1)*NpX; %index counter for vectorized vine structure
            
            %check if vines have just become latent, if so calculate mu_L 
            %and flip their LatentSwitch so that the calc only happens 1 time.
            if((vine(cnt).L(t-1) > 1e-8) && (vine(cnt).LatentSwitch == false))
                vine(cnt).mu_L = latentperiod(t,dt,Nsteps,mu_L_target,...
                    zeros(size(T)),T);
                vine(cnt).LatentSwitch = true;
            end

            % set initial conditions for time integration (could use deal here)
            y0(1) = vine(cnt).B(t-1); %(amount of population, surface area, that is berries)
            y0(2) = vine(cnt).P(t-1); %(total population, surface area, including berries and leaves)
            y0(3) = vine(cnt).S(t-1); %(initial susceptible population fraction)
            y0(4) = vine(cnt).L(t-1); %(initial latent population fraction)
            y0(5) = vine(cnt).I(t-1); %(initial infectious population fraction)
            y0(6) = vine(cnt).R(t-1); %(initial recovered population fraction)
            y0(7) = vine(cnt).E(t-1); %(initial external population fraction)
            y0(8) = vine(cnt).F(t-1); %(size of the spreading population)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
            DepFlux_sum_cnt = DepFlux_sum(cnt);
            vine_mu_L = vine(cnt).mu_L;

             [y] = TimeInt(odefun,t,dt,y0,DepFlux_sum_cnt,vine_mu_L);

            
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % set outputs
            vine(cnt).B(t) = y(1);
            vine(cnt).P(t) = y(2);
            vine(cnt).S(t) = y(3);
            vine(cnt).L(t) = y(4);
            vine(cnt).I(t) = y(5);
            vine(cnt).R(t) = y(6);
            vine(cnt).E(t) = y(7);
            vine(cnt).F(t) = y(8);

            %define a threshold for dispersal to start (equiv to 2.5mm diameter
            %sporulating colony (Calonnec et al)
            if(vine(cnt).I(t)>=0.25^2/4*pi/A)
                vine(cnt).IsInfect=true;
            end
            %turn off dispersal if we fall below the above size
            if(vine(cnt).I(t)<0.25^2/4*pi/A && vine(cnt).IsInfect==true)
                vine(cnt).IsInfect=false;
            end
        end
    end
    %% Scouting Routine
    if scoutingMethod == 3 && findSwitch ==0
        if t>=48 && mod(t,2)==0 % Starting on hour 48 we search with one drone every other hour
            speed = .04+1/600;  % Speed needed to search exactly 150 plants
            cost = cost + 100;  % each time one drone runs for one hour
            [infects,infectsFound] = Scouting(speed,vine,t,scoutingMethod,1,infects);
            if (t-1)/24 > 10 && mod(t,24)==0 % if we pass day 10 1000 dollars added per day
                cost = cost + 1000;
            end

        end
    else 
       if mod(t,24) == 0 
            switch scoutingMethod % speed and number of drones determined from Scouting_Opt for each method
                case 1 
                    amt = 4;
                    speed = .045;
                case 2
                    amt = 3;
                    speed = .055;
            end
            [infects,infectsFound] = Scouting(speed,vine,t,scoutingMethod,amt,infects);
            cost = cost + amt*100;
            if (t-1)/24 > 10
                cost = cost + 1000;
            end
       end
    end
    if infectsFound == 1 && findSwitch == 0 % the first time we find the infection we display it 
        tFound = t;
        disp('Infection Found')
        if onlyScout 
            return
        end
    end   

end

end

%% functions
function [y] = TimeInt(odefun,t,dt,y0,DepFlux_sum_cnt,vine_mu_L)
% trapizodal integration estimate
    y1 = odefun(t-1,y0,DepFlux_sum_cnt,vine_mu_L);
    y2 = odefun(t,y0+y1*dt,DepFlux_sum_cnt,vine_mu_L);
    y = y0 + ((y2+y1)/2)*(dt);           
end

function [infects,infectsFound] = Scouting(speed,vine,t,opts,amt,infects)
%inputs-Speed of drone,Vine Data from the simulation, current time, Pick scouting method | 1:Random search | 2: Modified random search | 3: grid search , number of drones, empty 50x50 matrix
%Outputs, 50x50 matrix of location of the infected plant, if the infection was located 
    global NpX NpY
    infectsFound = 0;
    DetectSize = (20*speed/10)^2/4*pi/5000;
    distMax = speed*3600;
    distUsed = 0;
    %% Random Search
    if opts == 1
        for a = 1:amt
            currLoc = [0,0]; % each drone starts in the corner 
            while distUsed < distMax && infectsFound ~= 1
                RandSearch = randi(NpX*NpY); % picking next plant
                distUsed = distUsed + sqrt((vine(RandSearch).X - currLoc(1))^2 + (vine(RandSearch).Y - currLoc(2))^2); % Calculating distance to next point and adding it to the total
                if distUsed > distMax % verifying target is reachable
                    break
                end
                if vine(RandSearch).I(t) >= DetectSize % checking the current plant
                    infects(vine(RandSearch).X+0.5,vine(RandSearch).Y+0.5) = 1;
                    infectsFound = 1;
                    return
                end
                currLoc = [vine(RandSearch).X,vine(RandSearch).Y]; % setting current location 

            end
        end
    end
    %% Modified Random Search
    if opts == 2
        for a = 1:amt
            corner = randi([1,4]);
            switch(corner)
                case 1
                    currLoc = [0,0];
                case 2
                    currLoc = [0,50];
                case 3 
                    currLoc = [50,0];
                case 4
                    currLoc = [50,50];
            end
            while distUsed < distMax && infectsFound ~= 1
                newrandx = randi([round(currLoc(1))-10, uptox(round(currLoc(1))+10)],"uint16");
                newrandy = randi([round(currLoc(2))-10, uptoy(round(currLoc(2))+10)],"uint16");
                while newrandy == 0 && newrandx == 0 
                    newrandx = randi([round(currLoc(1))-10, uptox(round(currLoc(1))+10)],"uint16");
                    newrandy = randi([round(currLoc(2))-10, uptoy(round(currLoc(2))+10)],"uint16");    
                end
                RandSearch = round(newrandx)+floor(newrandy)*50;
                distUsed = distUsed + sqrt((vine(RandSearch).X - currLoc(1))^2 + (vine(RandSearch).Y - currLoc(2))^2);

                if distUsed > distMax
                    break
                end
                if vine(RandSearch).I(t) >= DetectSize
                    infects(vine(RandSearch).X+0.5,vine(RandSearch).Y+0.5) = 1;
                    infectsFound = 1;
                    return
                end
                currLoc = [vine(RandSearch).X,vine(RandSearch).Y];
            end
        end
    end
    %% Grid Search
    if opts == 3
        t = (t-44)/2; 
        currLocInx = 1+(t-2)*3*50;
        distUsed = .5;
        adding = true; 
        numits =1;
        while distUsed<distMax
            if mod(numits,50)==1
                if adding
                    currLocInx = currLocInx+49;
                else
                    currLocInx = currLocInx+51;
                end
                adding = not(adding);
            end
            
            if vine(currLocInx).I(2*t+44) >= DetectSize
                infects(vine(currLocInx).X+0.5,vine(currLocInx).Y+0.5) = 1;
                infectsFound = 1;
                return
            end
            numits=numits+1;
            distUsed = distUsed+1;
            if adding
                currLocInx = 1+currLocInx;
            else
                currLocInx = currLocInx-1;
            end
        end
    end
end

function max = uptox(int)

if int>50
   max = 50;
else
   max = int;
end

end
function max = uptoy(int)

if int>=50
   max = 49;
else
   max = int;
end

end
