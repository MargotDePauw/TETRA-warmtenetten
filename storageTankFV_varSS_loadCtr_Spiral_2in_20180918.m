% this function calculates the thermal behavior of a STORAGE
% TANK WITHOUTH BOUYANT JET MIXING at in and outlet. 
% It is completely based on the type 4 storage of transys
% Only one inlet and one outlet are assumed. Connections are at upper and
% lower finite volume. No losses through top and bottom are assumed. This
% model is also equiped with a control logic system to generate a
% loading/not loading demand signal.
%
% Als je zowel een primaire én secundaire kant hebt (en dus eigenlijk 4 
% aansluitingen zou hebben), dien je het verschil van primair en secundair
% debiet te berekenen en als input te geven, dat komt bij de berekeningen 
% toch exact op het zelfde neer.
%
% De PDE die achter het model zit, is heeeeel erg vereenvoudigd erin 
% gestoken. In previous versions of this model (i.e. without _varStepSize),
% a constant step size is assumed. However, this might lead to enormous
% numerical noise (oscillations), overruling all physical meaning.
% Therfore, a minimum stepsize is applied in this model, that ensures a
% pre-set range of allowed oscillation.
%
% INPUTS
%%%%%%%%%
% inPhys
%   T_in [C]     temperature of the ingoing water
%   T_in_cir     temp of ingoing water from circulation pipes
%   T_in_spi [C]     temperature of the ingoing water of spiral
%   T_env [C]    temperature of environment
%   m_dot_top [kg/s] mass flow rate of the water at top of the tank. 
%                Positive if flow from top to bottom, 
%                negative if flow from bottom to top.
%   m_dot_bot [kg/s] mass flow rate of the water at bottom of the tank. 
%                Positive if flow from top to bottom, 
%                negative if flow from bottom to top.
%                difference between m_dot_top and m_dot_bot is circulation
%                water
%   m_dot_spi [kg/s] mass flow rate of the water through spiral. Positive if
%                flow from top to bottom, negative if flow from bottom to
%                top.
%   T_sp [°C]   set punt temperature of finite volume n_TT, same dim as
%               n_TT
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
%   UA_spi [W/K]overall heat tr coeff of internal spiral
%   i_cir       position of inlet of circulation pipe connection. between 1
%               and length(T_fv)
%
% OUTPUTS
%%%%%%%%%
%   T_top [C]   temperature of top connection (=T_in if m_dot_top>0)
%   T_bot [C]   temperature of bottom connection (=T_in if m_dot_bot<0)
%   T_fv [C]    temperatures of finite volumes
%   T_spi_out [C] temp of water going out spiral
%   Q_dot_loss [W] total heat loss to surroundings
%   Q_dot_spi [W] total heat load recieved from spiral (negative if spiral
%                is being heated)
%   e_loading [BOOL] 1 if the vessel wants to be loaded
%
% freek.vanriet@uantwerpen.be

% adjustments
% 2017 07 12   error fixed: Q_dot_loss extracted from FVs, not added
% 2017 07 12   _varStepSize: doing subsimulations
% 2017 08 28   changing definition of T_in and T_out to T_top and T_bot
% 2017 09 01    loading control added (i.e. signal e_loading to production
%               system to tell whether or not loading is required)
% 2017 09 06    possibility to have two sensors added
% 2017 12 04    give always non-NaN T_top and T_bot
% 2018 09 12    add spiral
function outVars=storageTankFV_varSS_loadCtr_Spiral_2in_20180918(inPhys,inHis,param,timestep)
%% assign variables
% inPhys
%   see called function

% inHis
T=inHis.T_fv; %states of the finite volume of the storage tank
e_loading_old=inHis.e_loading;

% param
D=param.D;
L=param.L;
dT_max=param.dT_max;
n_TT=param.n_TT;

n_fv=length(T);


%% preliminary calculations
outVars=storageTankFiniteVolumes_Spiral_2in_20180918(inPhys,inHis,param,timestep);

%% estimating max temperature change and setting timestep
Q_dot_max=max(abs(outVars.Q_dot));
A=pi*(D/2)^2;%area of section between two finite volumes
C=4187*1000*L/n_fv*A;
dt_max=abs(dT_max*C/Q_dot_max);


%% simultions in sub time 

if timestep>dt_max % max temperature change will be exceeded --> need for 'subsimulations'
    n_subsim=ceil(timestep/dt_max);
    timestep_subsim=timestep/n_subsim;
    
   % initialise storage array for T_out and Q_dot ...
   T_top=zeros(n_subsim,1);
   T_bot=T_top;
   Q_dot_bulk=T_top;
   Q_dot_con=T_top;
   Q_dot_loss=T_top;
   Q_dot_spi=T_top;
   
   T_spi_out=zeros(n_subsim,n_fv);
   Q_dot=T_spi_out;

%    if n_subsim>100
%        keyboard
%    end
    for i=1:n_subsim
        %call
        outVars_subsim=storageTankFiniteVolumes_Spiral_2in_20180918(inPhys,inHis,param,...
            timestep_subsim);
        
        %inputs for next iteration
        %   inHis
        inHis.T_fv=outVars_subsim.T_fv;
        %   inPhys stays the same
        
        %outputs        
        T_top(i)=outVars_subsim.T_top;
        T_bot(i)=outVars_subsim.T_bot;
        T_spi_out(i,:)=outVars_subsim.T_spi_out;
        Q_dot_bulk(i)=outVars_subsim.Q_dot_bulk;
        Q_dot_con(i)=outVars_subsim.Q_dot_con;
        Q_dot_loss(i)=outVars_subsim.Q_dot_loss;
        Q_dot_spi(i)=outVars_subsim.Q_dot_spi;
        Q_dot(i,:)=outVars_subsim.Q_dot;
%         keyboard
    end
    outVars.Q_dot_bulk=mean(Q_dot_bulk);
    outVars.Q_dot_con=mean(Q_dot_con);
    outVars.Q_dot_loss=mean(Q_dot_loss);
    outVars.Q_dot_spi=mean(Q_dot_spi);
    outVars.Q_dot=mean(Q_dot,1);
    outVars.T_fv=outVars_subsim.T_fv;
    outVars.T_top=mean(T_top);
    outVars.T_bot=mean(T_bot);
    outVars.T_spi_out=mean(T_spi_out,1);
    
%     % assign output temperatures
%     if m_dot>0 %flow from top to bottom/ loading
%         outVars.T_bot=T_out_mean;
%     elseif m_dot<0 %flow from bottom to top/deloading
%         outVars.T_top=T_out_mean;
%     else %no flow
%         outVars.T_top=NaN;
%         outVars.T_bot=NaN;
%     end

end

%% loading/not loading control signal
if length(n_TT)==1    
    T_sp=inPhys.T_sp;
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
    T_sp_top=inPhys.T_sp(n_TT==n_TT_top);
    T_sp_bot=inPhys.T_sp(n_TT==n_TT_bot);
    
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
    displ('I do not have enough probes to put all your damn sensors in me, you idiot! Please reduce the number of sensors to a maximum of two')
    error('Abort mission')
end

outVars.e_loading=e_loading;
end