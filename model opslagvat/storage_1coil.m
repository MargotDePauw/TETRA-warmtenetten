% this function calculates the thermal behavior of a STORAGE
% TANK with 1 or 2 coils , based on the sat boiler model used in install2020 
% up to 2 sensors possible
% flow in HE is input parameter
%
% INPUTS
%%%%%%%%%
% inPhys
%   T_in [C]        temperature of the ingoing water
%   T_env [C]       temperature of environment
%   m_dot [kg/s]    mass flow rate of the water through the tank. Positive if
%                   flow from top to bottom, negative if flow from bottom to
%                   top.
%   T_sp [°C]       set punt temperature of finite volume n_TT, same dim as
%                   n_TT
%   T_in_HE [C]     ingoing temperature heat exchanger
%   m_dot_HE [kg/s] mass flow heat exchanger
%
% inHis
%   T_fv [C]    historical temperatures of all finite volumes
%   e_loading   1 if the vessel wanted to be loaded the previous moment in
%               time
%
% param
%   D [m]       inner diameter of tank
%   L [m]       height of the tank
%   UA [W/K]    overall heat transfer coefficient
%   dT_max [°C] maximum absolute change in temperature 
%               (ensures little oscillation)
%   n_TT        at partial volume(s) n_TT, the temperature sensor(s) is/are
%               located. Needed to produce a control signal for the loading
%               1x1 or 1x2
%   DT_hys [°C] Delta T for hysteresis control, if dimensions of n_TT equal
%               to 1x2, DT_hys doesn't matter
%   R1 [K/W]     R-value of coil (1 loop)
%  coil heat exchanger is supposed to have 12 loops of 1m length and
%  filling bottom half of storage vessel, so if we have 50 voluminas, 
%  each volumina gets heat from circa 0.5m loop
%
% OUTPUTS
%%%%%%%%%
%   T_top [C]   temperature of top connection (=T_in if m_dot>0)
%   T_bot [C]   temperature of bottom connection (=T_in if m_dot<0)
%   T_fv [C]    temperatures of finite volumes
%   Q_dot_loss [W] total heat loss to surroundings
%   e_loading [BOOL] 1 if the vessel wants to be loaded
%   T_ret_HE [C]    outgoing temperature heat exchanger
%
% freek.vanriet@uantwerpen.be

% adjustments
% 2017 07 12   error fixed: Q_dot_loss extracted from FVs, not added
% 2017 07 12   _varStepSize: doing subsimulations
% 2017 08 28   changing definition of T_in and T_out to T_top and T_bot
% 2017 09 01    loading control added (i.e. signal e_loading to production
%               system to tell whether or not loading is required)
% 2017 09 06    possibility to have two sensors added
% 2017 12 13    heat exchanger added to load storage vessel

function outVars=sat_boiler_20180115(inPhys,inHis,param,timestep)
%% assign variables
% inPhys
T_in=inPhys.T_in; %supply tank
T_env=inPhys.T_env;
m_dot=inPhys.m_dot;%water through tank
T_in_HE1 = inPhys.T_in_HE; % temp of fluid entering heat exchanger1

m_dot_HE = inPhys.m_dot_HE;%flow in HE1, calculated in previous timestep to use in this timestep


% inHis
T=inHis.T_fv; %states of the finite volume of the storage tank, starting from top layer
e_loading_old=inHis.e_loading;

% param
D=param.D;
L=param.L;
UA=param.UA;%W/K
dT_max=param.dT_max;
n_TT=param.n_TT;
R_HE = param.R_HE;
n_coil_start=min(param.n_HE);%layer where HE enters storage vessel
n_coil_end=max(param.n_HE);%layer where HE leaves storage vessel
n_fv=length(T);



%% preliminary calculations of incoming temperatures bulk flow for each layer
if m_dot>0 %flow from top to bottom/ loading
    disp('loading storage by direct flow not possible in this model');
    T_bulk_in=[T_in, T(1:end-1)];
    outVars.T_top=T_in;
    outVars.T_bot=T(end);
elseif m_dot<0 %flow from bottom to top/deloading
    T_bulk_in=[T(2:end),T_in];
    outVars.T_top=T(1);
    outVars.T_bot=T_in;
else %no flow
    T_bulk_in=T; % doesn't matter, because it will be multiplied by zero
    outVars.T_top=T(1);
    outVars.T_bot=T(end);
end

%% preliminary calculations conduction between layers
%k=0.65; %W/mK thermal conductivity of water at 60°C
A=pi*(D/2)^2;%area of section between two finite volumes
delta_k=50*(pi*D*0.02)/A;%conduction through tank wall
k=0.65;%+delta_k;
T_cond_up=[T(1),T(1:end-1)];%no conduction at top
T_cond_low=[T(2:end),T(end)];%no conduction at bottom

%% heat losses to surroundings for each layer
Q_dot_loss=UA/n_fv*(T-T_env);%heat loss to surroundings

%% heat flow from coil

Q_HE = zeros(1,n_fv);
T_coil = zeros(1,n_coil_end-n_coil_start+2);%entering value for each layer with coil and outgoing value of last layer zeros(1,n_fv-n_coil_start+2)
if e_loading_old == 1 %flow in coil
T_coil(1,1)=T_in_HE1;
for i= 1 : n_coil_end-n_coil_start+1  %n_fv-n_coil_start+1
    %R_HE = 1/UA for 1m loop (=length in each layer) (*2 if only 0.5m, etc)
    T_coil(1,i+1) = T_coil(1,i)+(T(n_coil_start+i-1)-T_coil(1,i))*(1-exp(-1/(R_HE*m_dot_HE*4187)));
end
Q_HE(n_coil_start:n_coil_end) = -diff(T_coil)*transpose(m_dot_HE)*4187; %Q_HE(n_coil_start:n_fv)
end

%% energy balance for each layer
Q_dot=4187*abs(m_dot)*(T_bulk_in-T)... %energy flow due to bulk flow
      + k*A/(L/n_fv)*(T_cond_up-2*T+T_cond_low)... %conduction
      - Q_dot_loss ... % losses to surroundings
    + Q_HE;%from coin
%% estimating max temperature change and setting timestep
Q_dot_max=max(abs(Q_dot));
C=4187*1000*L/n_fv*A;
dt_max=abs(dT_max*C/Q_dot_max);


%% simultions in sub time 

if timestep<=dt_max % max temperature change OK
    outVars.Q_dot_loss=sum(Q_dot_loss);
    outVars.T_fv=T+Q_dot*timestep/(4187*1000*L/n_fv*A);
    outVars.Toutcoil = T_coil(end);%0 when no flow in HE
else % max temperature change will be exceeded --> need for 'subsimulations'
    disp('subsimulations');
    n_subsim=ceil(timestep/dt_max);
    timestep_subsim=timestep/n_subsim;
    
   % initialise storage array for T_out and Q_dot_loss
   
   T_out=zeros(n_subsim,1);
   Q_dot_loss=T_out;
   Toutcoil = T_out;
   
    for i=1:n_subsim
        %call
        outVars_subsim=storageTankFiniteVolumes_20171218(inPhys,inHis,param,...
            timestep_subsim);
        
        %inputs for next iteration
        %   inHis
        inHis.T_fv=outVars_subsim.T_fv;
        %   inPhys stays the same
        
        %outputs        
        T_out(i)=outVars_subsim.T;
        Q_dot_loss(i)=outVars_subsim.Q_dot_loss;
        Toutcoil(i) = outVars_subsim.Toutcoil;
    end
    outVars.Q_dot_loss=mean(Q_dot_loss);
    outVars.T_fv=outVars_subsim.T_fv;
    T_out_mean=mean(T_out);
    outVars.Toutcoil = mean(Toutcoil);
    
    % assign output temperatures
    if m_dot>0 %flow from top to bottom/ loading
        outVars.T_top=outVars_subsim.T_fv(1);
        outVars.T_bot=T_out_mean;
    elseif m_dot<0 %flow from bottom to top/deloading
        outVars.T_top=T_out_mean;
        outVars.T_bot=outVars_subsim.T_fv(end);
    else %no flow
        outVars.T_top=outVars_subsim.T_fv(1);
        outVars.T_bot=outVars_subsim.T_fv(end);
    end
end


%%reversion elimination algorithm !!!!!!!!evt terug weglaten
Ttemp = outVars.T_fv;
ind1 = find(diff(Ttemp)>0.001);%indices where next (lower) temp is higher
while not(isempty(ind1))
ind2 = ind1+1;
T1=Ttemp(ind1);
T2=Ttemp(ind2);
T12=cat(1,T1,T2);
Tgem= mean(T12,1);
rep0=[1:length(Tgem)];
rep=sort([rep0,rep0]);
Tgemd=Tgem(rep);
ind = sort([ind1,ind2]);
Ttemp(ind)= Tgemd;
%check energy
if (mean(Ttemp) - mean(outVars.T_fv))<0.00001
outVars.T_fv=Ttemp;
else
    disp('fout in reversion elimination algorithm')
end
ind1 = find(diff(Ttemp)>0.001);
end




%% loading/not loading control signal
if length(n_TT)==1    
    T_sp=param.T_sp;%inPhys.T_sp;
    DT_hys=param.DT_hys;
    switch e_loading_old
        case 0 % previous time the vessel did not have any desires to be warmed
               % up (he felt happy)
            if outVars.T_fv(n_TT)<T_sp-DT_hys;
                e_loading=1; % not so happy --> need warmth
            else
                e_loading=0; % still happy
            end
        case 1 % previous time, the vessel was warming up. Is he happy now?
            if outVars.T_fv(n_TT)>T_sp+DT_hys
                e_loading=0; % NO! It's to hot, I cannot take this!
            else
                e_loading=1; % I don't no... Give me some extra warmth please.
            end
    end
    
elseif length(n_TT)==2
    n_TT_top=min(n_TT);
    n_TT_bot=max(n_TT);
    T_sp_top=param.T_sp(n_TT==n_TT_top);%inPhys.T_sp(n_TT==n_TT_top);
    T_sp_bot=param.T_sp(n_TT==n_TT_bot);
    
    switch e_loading_old
        case 0 % previous time the vessel did not have any desires to be warmed
               % up (he felt happy)
            if outVars.T_fv(n_TT_top)<T_sp_top && outVars.T_fv(n_TT_bot)<T_sp_bot
                e_loading=1; % not so happy --> need warmth
            else
                e_loading=0; % still happy
            end
        case 1 % previous time, the vessel was warming up. Is he happy now?
            if outVars.T_fv(n_TT_top)>T_sp_top && outVars.T_fv(n_TT_bot)>T_sp_bot
                e_loading=0; % NO! It's to hot, I cannot take this!
            else
                e_loading=1; % I don't no... Give me some extra warmth please.
            end
    end

else
    disp('I do not have enough probes to put all your damn sensors in me, you idiot! Please reduce the number of sensors to a maximum of two')
    error('Abort mission')
end

outVars.e_loading=e_loading;
outVars.T_ret = outVars.Toutcoil; % Toutcoil =0 if no flow in HE...

outVars.T_out = outVars.T_top;
outVars.T_in = inPhys.T_in;

outVars.Q_loss = -outVars.Q_dot_loss;%W
outVars.P_HE = 4.185*m_dot_HE*(inPhys.T_in_HE-outVars.T_ret);%power to storage
outVars.P_out = 4.185*m_dot*(outVars.T_out-outVars.T_in);
outVars.T_fv_mean = mean(outVars.T_fv);

end