% DESCRIPTION: 
% This function simulates the dynamic thermal behaviour of a solar thermal
% collector. The model is based on:
% TRNSYS Type 832 v5.00 �Dynamic Collector Model by Bengt Perers�
% with no wind speed and IR dependant losses.
% The internal energy change is calculated with the mean temperature of in-
% and outlet.
%
% INPUTS: 
%   'inPhys': physical inputs
%           m_dot(kg/s)     mass flow through collector
%           T_in(�C)        fluid temperature at inlet
%           T_amb(�C)       air temperature outdoor
%           theta (�)       angle of incidence
%           azimuth_sol (�)
%           zenith_sol (�)
%           I_b (W/m�)      beam radiation on the collector
%           I_d (W/m�)      diffuse radiation (gnd + sky) on the collector
%
%   'inHis': historical values (for, among other reasons,
%                 solving differential equations)
%           T_m(�C)       mean temp. at outlet of previous timestep ((T_in+T_out)/2)
%           T_out(�C)       outlet temp at previous timestep
%           
%
%   'param': parameters
%           C(J/�C)          overall thermal capacity (i.e. capacitance)
%           A(m�)            aperture area of collector
%           k0(-)            zero-loss efficientie (FR*(tau*alfa)n)
%           k1(W/m�K)        linear heat loss coeff
%           k2(W/m�K�)       quadratic heat loss coeff
%           c(J/kg�C)        specific heat capacity of collector fluid
%           m_dot_test (kg/s)    flow rate at test
%           beta_collect (degrees)   angle between collector and horizontal (slope)
%           IAM              Incidence Angle Modifier data structure of
%                            beam irradiance. 'l' and 't' are different for
%                            e.g. vacuum tubes. If only a single (set of)
%                            IAM available (e.g. flat plate) just take
%                            either 'l' or 't' for unity.
%               .theta_l     longitudinal incidence angle (0�... 90�)
%               .IAM_l       longitudinal IAM
%               .theta_t     transversal incidence angle (0�... 90�)
%               .IAM_t       transversal IAM
%           azimuth_collect
%           T_out_max(�C)    maximum out temperature -> shutdown
%           T_in_max(�C)     maximum in temperature -> shutdown
%           ro_g (-)         ground reflectance
%
%
%   1x1 'timestep'
%           timestep(s)      simulation timestep
%
% OUTPUTS:
%   'outVars': outputs of model
%           T_out(�C)       fluid temp. at outlet of new timestep. This 
%                           should not be used.
%           T_out_mean(�C)  mean outlet temperature during timestep. This
%                           should be used as input for next component
%           T_m             arithmic mean fluid temperature of new
%                           timestep. This should be used as initial
%                           condition of next timestep.
%           T_m_mean        mean of arithmic mean fluid temperature during 
%                           timestep. This shoul be used for energy balance
%                           checks and within this function it is used for
%                           heat loss calculations (Q_dot_loss)
%           Q_dot_loss(W)   mean heat flux to surroundings
%
% freek.vanriet@uantwerpen.be; 19/03/2020
% based on boiler_tempControl_20170926
% 
% adjustments
%28/05/2020 add calculation of global IAM 
% only 1 parameter structure IAM (instead of separate ones for beam and
% diffuse)
% add parameter beta (slope of collector)and azimuth of collector (for
% calculation of thata l and t)
% add solar azimuth and zenith as input (for calculation of theta l and t)
% add parameter ground reflectance for calculation of global IAM
% add parameter k0, A
% inHis: add T_out, T_run
% diffuse radiation sky and ground as separate inputs
%add correction for flow rate different from test (so mtest as new
%parameter)


function outVars=solCol_20200318(inPhys,inHis,param,timestep)
%% assign function inputs

%inPhys
     m_dot=inPhys.m_dot;
     T_in=inPhys.T_in;
     T_amb=inPhys.T_amb;
     I_b=inPhys.I_b;
     I_dsky=inPhys.I_dsky;
     I_dgnd=inPhys.I_dgnd;
     azimuth_s=inPhys.azimuth_sol;
     zenith_s=inPhys.zenith_sol;
     theta=inPhys.theta;
    
%inHist
    T_out=inHis.T_out;
   
    
%param
    C=param.C;
    A=param.A;
    k0=param.k0;
    k1=param.k1;
    k2=param.k2;
    c=param.c_water;
    m_dot_test=param.m_dot_test;
    azimuth_c=param.azimuth_collect;
    beta = param.beta_collect;
    ro_g=param.ro_g;
    IAM.theta_t=param.IAM.theta_t;
    IAM.IAM_t=param.IAM.IAM_t;
    IAM.theta_l=param.IAM.theta_l;
    IAM.IAM_l=param.IAM.IAM_l;

%% control

% not applicable within this model?
%% IAMs and radiative balance

% projections of incidence angle on longitudinal and transversal plane
%           theta_l (�)     longitudinal incidance angle of sunlight on
%                           collector plane
%           theta_t (�)     transversal incidance angle of sunlight on
%                           collector plane
if zenith_s == 90 
    theta_l=90;
    theta_t=90;
else
    if theta==90
        theta_l=90;
        theta_t=90;
    else
        theta_l= abs(beta-atand(tand(zenith_s)*cosd(azimuth_c-azimuth_s)));
        theta_t= abs(atand((sind(zenith_s)*sind(azimuth_c-azimuth_s))/(cosd(theta))));%abs because of neg angles, is it ok?
    end
end
    

% IAMs
%   prevent extrapolation
theta_l= min(max(theta_l, IAM.theta_l(1)),IAM.theta_l(end));
theta_t= min(max(theta_t, IAM.theta_t(1)),IAM.theta_t(end));



%   interpolate

%IAM for beam radiation
IAMb=interpn(IAM.theta_t,IAM.IAM_t,theta_t)*interpn(IAM.theta_l,IAM.IAM_l,theta_l);
%IAM for diffuse radiation
%equivalent incident angle for diffuse radiation
theta_sky = 59.68-0.1388*beta+0.001497*beta*beta;
theta_gnd = 90-0.5788*beta+0.002693*beta*beta;
IAMdsky = interpn(IAM.theta_t,IAM.IAM_t,theta_sky)*interpn(IAM.theta_l,IAM.IAM_l,theta_sky);  %is it correct to take product here,?
IAMdgnd = interpn(IAM.theta_t,IAM.IAM_t,theta_gnd)*interpn(IAM.theta_l,IAM.IAM_l,theta_gnd);  %is it correct to take product here,?
%global IAM
I_t = I_b+I_dsky+I_dgnd;
if I_t > 0
%    IAM = (I_b*IAMb+(1+cosd(beta))*0.5*I_t*IAMdsky+ro_g*(1-cosd(beta))*0.5*I_t*IAMdgnd)/I_t;  %or is it better to get sky diffuse and ground diffuse radiation as input??
    IAM = (I_b*IAMb+I_dsky*IAMdsky+I_dgnd*IAMdgnd)/I_t;
else
    IAM=1;
end
% 

    

if C == 0  %steady state model
    if m_dot == 0
        m_dot_temp = m_dot_test; %enables calculation of T_out as control value
    else
        m_dot_temp=m_dot;
    end
    %correction for flow rate different from test (only for steady state calculation)
    f= - m_dot_test*log(1-(k1*A/(m_dot_test*c)));
    corr_flow = (m_dot_temp/m_dot_test)*(1-exp(-f/m_dot_temp))/(1-exp(-f/m_dot_test));
    DT=(T_in+T_out)*0.5-T_amb;  %efficiency (heat losses) is based on T_out from previous timestep
    if I_t >0
        eff = (k0*IAM-k1*DT/I_t-k2*DT*DT/I_t)*corr_flow;
        T_out_new= T_in + eff*I_t*A/(c*m_dot_temp); %constant value in this timestep, cf no capacity
    else
        T_out_new=T_in + (-k1*DT*A-k2*DT*DT*A)*corr_flow/(c*m_dot_temp); %or Nan??
    end
    T_out_mean=T_out_new;
    Q_dot=c*m_dot*(T_out_new-T_in);
    %DT=(T_in + T_out_mean)*0.5-T_amb;  %new DT with new Tout
    Q_dot_loss_new=(k1*DT+k2*DT*DT)*corr_flow*A;  
    Q_dot_rad = k0*IAM*A*I_t;
end



%% dynamics

% ADJUST IF NECESSARY: DEFINITION OF MEAN TEMEPRATURE
if C > 0  %model with collector capacitance
    % test coefficients k0,k1,k2 should be corrected here, as the
    % differential equation is with F' instead of
    % FR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % radiative balance
    Q_dot_rad = k0*IAM*A*I_t; 
    k1_corr = k1*A;
    k2_corr=k2*A;
    % make it possible to have T_in being NaN as input if no flow
    if m_dot<=0 % the model should be able to have T_in as a NaN
        T_in=0;    % NO PHYSICAL MEANING
    end

    DT_0=(T_in+T_out)/2-T_amb;%DT_0=(T_in+T_out_0)/2-T_amb;

    if DT_0>5 % OK to use squared term
        % define parameters of "diff(DT(t), t) = -a2*DT(t)^2 - a1*DT(t) + a0"
        a2=k2_corr/C;
        a1=(k1_corr+2*c*m_dot)/C;
        a0=(Q_dot_rad-2*c*m_dot*(T_amb-T_in))/C;%+ replaced with -2*c*m_dot enz
        
        DT = real((tanh(timestep*sqrt(4*a2*a0 + a1^2)/2 + atanh((2*DT_0*a2 + a1)/sqrt(4*a2*a0 + a1^2)))*sqrt(4*a2*a0 + a1^2) - a1)/(2*a2)); %real, because result of atanh has complex NOTATION
        
    else % NOT OK: DT can switch between + and -
        % define parameters of "diff(DT(t), t) = - a1*DT(t) + a0"
        % and fit linear curve through sec. order polynomial
        DT_fit = DT_0+0.000001; %????????, is DT_0 ok as value for DT_fit (0.00001 otherwise Nan is DT_0=0)
        Q_dot_loss_fit=k1_corr*DT_fit+sign(DT_fit)*k2_corr*DT_fit^2;
        k_=Q_dot_loss_fit/DT_fit;%was DT_fit/Q_dot_loss_fit; ????????????????????????????????????????????????????????????????????????
        a1_=(k_+2*c*m_dot)/C;
        a0=(Q_dot_rad-2*c*m_dot*(T_amb-T_in))/C;%+ replaced with -2*c*m_dot enz
        
        DT = a0/a1_ + exp(-a1_*timestep)*(DT_0 - a0/a1_);
       % a1=a1_; %for calculation T_out_mean
    end
    T_out_new =(DT+T_amb)*2-T_in;
    %T_out_mean= (T_out + T_out_new)*0.5;  %NOT 100% CORRECT!!! (T_out is not lineair function in time)
    %approximation by exponential curve instead of lineair, 
    %exponential curve is solution of diff (DT,t)= -b1DT +b0 (same lineair approx as for DT_0<5)
    b0=a0;
    b1= ((k1_corr*DT_0+sign(DT_0)*k2_corr*DT_0^2)/(DT_0+0.0000001)+2*c*m_dot)/C;
    %T_out_mean = (b0*timestep-T_out_new+T_out)/(b1*timestep); %approximation by exponential curve instead of linair
    DT_mean = (b0*timestep-DT+DT_0)/(b1*timestep);
    T_out_mean=(DT_mean+T_amb)*2-T_in;
    
    Q_dot=c*m_dot*(T_out_mean-T_in);
    %DT_mean = (T_in+T_out_mean)/2-T_amb;
    if DT_0>5
        Q_dot_loss_new=(k1*DT_mean+k2*DT_mean*abs(DT_mean))*A; 
    else
       Q_dot_loss_new =k_*DT_mean;
    end
    corr_flow=1; % no correction needed here
end
%% calculate other outputs and assign outputs to output structure


% assign outputs
outVars.T_out=T_out_new;
outVars.T_out_mean=T_out_mean;
outVars.Q_dot_loss=Q_dot_loss_new;%in W
outVars.Q_dot=Q_dot;
outVars.Q_dot_rad=Q_dot_rad;
outVars.IAM=IAM;
outVars.IAMb=IAMb;
outVars.corr_flow=corr_flow;
outVars.m_dot=m_dot;
outVars.Tcoll=DT+T_amb; % collector temp (for energy balance)
outVars.DT_0=DT_0; %for testing of energy balance, difference in case of lineair and kwadratic losses
end