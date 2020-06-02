% DESCRIPTION: 
% This function simulates the dynamic thermal behaviour of a solar thermal
% collector. The model is based on:
% TRNSYS Type 832 v5.00 „Dynamic Collector Model by Bengt Perers“
% with no wind speed and IR dependant losses.
% The internal energy change is calculated with the mean temperature of in-
% and outlet.
%
% INPUTS: 
%   'inPhys': physical inputs
%           m_dot(kg/s)     mass flow through collector
%           T_in(°C)        fluid temperature at inlet
%           T_amb(°C)       air temperature outdoor
%           theta (°)       angle of incidence
%           azimuth zon?
%
%   'inHis': historical values (for, among other reasons,
%                 solving differential equations)
%           T_m(°C)       mean temp. at outlet of previous timestep ((T_in+T_out)/2)
%
%   'param': parameters
%           C(J/°C)          overall thermal capacity (i.e. capacitance)
%           k1(W/m²K)        linear heat loss coeff
%           k2(W/m²K²)       quadratic heat loss coeff
%           c(J/kg°C)        specific heat capacity of collector fluid
%           IAMb             Incidence Angle Modifier data structure of
%                            beam irradiance. 'l' and 't' are different for
%                            e.g. vacuum tubes. If only a single (set of)
%                            IAM available (e.g. flat plate) just take
%                            either 'l' or 't' for unity.
%               .theta_l     longitudinal incidence angle (0°... 90°)
%               .IAM_l       longitudinal IAM
%               .theta_t     transversal incidence angle (0°... 90°)
%               .IAM_t       transversal IAM
%           IAMd             Incidence Angle Modifier data structure of
%                            diffuse irradiance. Typically both IAM_l and
%                            IAM_t equal to unity??
%               .theta_l     longitudinal incidence angle (0°... 90°)
%               .IAM_l       longitudinal IAM
%               .theta_t     transversal incidence angle (0°... 90°)
%               .IAM_t       transversal IAM
%           azimuth collector?
%           T_out_max(°C)    maximum out temperature -> shutdown
%           T_in_max(°C)     maximum in temperature -> shutdown
%
%   1x1 'timestep'
%           timestep(s)      simulation timestep
%
% OUTPUTS:
%   'outVars': outputs of model
%           T_out(°C)       fluid temp. at outlet of new timestep. This 
%                           should not be used.
%           T_out_mean(°C)  mean outlet temperature during timestep. This
%                           should be used as input for next component
%           T_m             arithmic mean fluid temperature of new
%                           timestep. This should be used as initial
%                           condition of next timestep.
%           T_m_mean        mean of arithmic mean fluid temperature during 
%                           timestep. This shoul be used for energy balance
%                           checks and within this function it is used for
%                           heat loss calculations (Q_dot_loss)
%           Q_dot_loss(W)   mean heat flux to surroundings
%           t_run (s)       new calculated run/off time
%
% freek.vanriet@uantwerpen.be; 19/03/2020
% based on boiler_tempControl_20170926
% 
% adjustments
%

function outVars=solCol_20200318(inPhys,inHis,param,timestep)
%% assign function inputs

%inPhys
     m_dot=inPhys.m_dot;
     T_in=inPhys.T_in;
     T_amb=inPhys.T_amb;
     I_b=inPhys.I_b;
     I_d=inPhys.I_d;
    
%inHist
    T_out=inHis.T_out;
    
%param
    C=param.C;
    k1=param.k1;
    k2=param.k2;
    c=param.c_water;

%% control

% not applicable within this model?
%% IAMs and radiative balance

% projections of incidence angle on longitudinal and transversal plane
%           theta_l (°)     longitudinal incidance angle of sunlight on
%                           collector plane
%           theta_t (°)     transversal incidance angle of sunlight on
%                           collector plane

% IAMs
%   prevent extrapolation

%   interpolate

% 

% radiative balance

%   Q_dot_rad = ...        

%% dynamics    

% ADJUST IF NECESSARY: DEFINITION OF MEAN TEMEPRATURE

% make it possible to have T_in being NaN as input if no flow
if m_dot<=0 % the model should be able to have T_in as a NaN
    T_in=0;    % NO PHYSICAL MEANING
end

DT_0=(T_in+T_out_0)/2-T_amb;

if DT_0>5 % OK to use squared term
    % define parameters of "diff(DT(t), t) = -a2*DT(t)^2 - a1*DT(t) + a0"
    a2=k2/C;
    a1=(k1+2*c*m_dot)/C;
    a0=(Q_dot_rad+2*c*m_dot*(T_amb-T_in))/C;
    
    DT = (tanh(timestep*sqrt(4*a2*a0 + a1^2)/2 + arctanh((2*DT_0*a2 + a1)/sqrt(4*a2*a0 + a1^2)))*sqrt(4*a2*a0 + a1^2) - a1)/(2*a2);
    
else % NOT OK: DT can switch between + and -
    % define parameters of "diff(DT(t), t) = - a1*DT(t) + a0"
    % and fit linear curve through sec. order polynomial
    Q_dot_loss_fit=k1*DT_fit+sign(DT_fit)*k2*DT_fit^2;
    k_=DT_fit/Q_dot_loss_fit;
    a1_=(k_+2*c*m_dot)/C;
    a0=(Q_dot_rad+2*c*m_dot*(T_amb-T_in))/C;
    
    DT = a0/a1_ + exp(-a1_*timestep)*(DT_0 - a0/a1_);
    
end

%% calculate other outputs and assign outputs to output structure


% time running or off 
if Q_dot_con==0 % HP is off
    if t_run<0 % HP was off
        t_run_new=t_run-timestep; % is off for a longer time
    else % HP was on
        t_run_new=-timestep; %is off since one timestep
    end
else % HP is on
    if t_run<0 %  was off
        t_run_new=timestep; % is on since one timestep
    else %  was on
        t_run_new=t_run+timestep; % is on for a longer time
    end
end

% assign outputs
outVars.T_out=T_out_new;
outVars.T_out_mean=T_out_mean;
outVars.Q_dot_loss=Q_dot_loss_new;%in W
outVars.t_run=t_run_new;
outVars.Q_dot_con=Q_dot_con;
outVars.A_dot=A_dot;

end