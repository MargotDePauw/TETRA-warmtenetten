% this function calculates the thermal behavior of a STORAGE
% TANK WITHOUTH BOUYANT JET MIXING at in and outlet. 
% It is completely based on thetype 4 storage of transys
% Only one inlet and one outlet are assumed. Connections or at upper and
% lower finite volume. No losses through top and bottom are assumed.
%
% Als je zowel een primaire én secundaire kant hebt (en dus eigenlijk 4 
% aansluitingen zou hebben), dien je het verschil van primair en secundair
% debiet te berekenen en als input te geven, dat komt bij de berekeningen 
% toch exact op het zelfde neer.
%
% De PDE die achter het model zit, is heeeeel erg vereenvoudigd erin 
% gestoken. Let dus op met te grote tijdstappen (wat is te groot? dat weet
% ik niet) en te kleine ‘finite volumes’ (wat is te klein? Dat weet ik ook 
% niet) . Dit laatste lijkt niet logisch, maar beide hebben te maken met 
% dat hoe groter de temperatuursverandering binnen 1 tijdstap is, hoe meer 
% mijn  berekening afwijkt van de eigenlijke PDE. Te grote finite volumes 
% leiden dan weer tot te veel uitsmering en te weinig stratificatie. 
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
%
% inHis
%   T_fv [C]    historical temperatures of all finite volumes
%
% param
%   D [m]       inner diameter of tank
%   L [m]       height of the tank
%   UA [W/K]    overall heat transfer coefficient
%   UA_spi [W/K]overall heat tr coeff of internal spiral
%   i_cir       position of inlet of circulation pipe connection. between 1
%               and length(T_fv)
%
% OUTPUTS
%%%%%%%%%
%   T_top [C]   temperature of top connection (=T_in if m_dot>0)
%   T_bot [C]   temperature of bottom connection (=T_in if m_dot<0)%   T_fv [C]    temperatures of finite volumes
%   T_spi_out [C] temp of water going out spiral
%   Q_dot_loss [W] total heat loss to surroundings
%   Q_dot_spi [W] total heat load recieved from spiral (negative if spiral
%                is being heated)
%
% freek.vanriet@uantwerpen.be

% adjustments
% 2017 07 12   error fixed: Q_dot_loss extracted from FVs, not added
% 2018 09 12    add spiral
% 2018 09 18    add extra inlet for circulation pipe of DHW

function outVars=storageTankFiniteVolumes_Spiral_2in_20180918(inPhys,inHis,param,timestep)
%% assign variables
% inPhys
T_in=inPhys.T_in;
T_in_cir=inPhys.T_in_cir;
T_in_spi=inPhys.T_in_spi;
T_env=inPhys.T_env;
m_dot_top=inPhys.m_dot_top;
m_dot_bot=inPhys.m_dot_bot;
m_dot_spi=inPhys.m_dot_spi;

% inHis
T=inHis.T_fv; %states of the finite volume of the storage tank

% param
D=param.D;
L=param.L;
UA=param.UA;
UA_spi=param.UA_spi;
i_cir=param.i_cir;

n_fv=length(T);

%% warning for timestep
if timestep>30 % depends actually on Q_dot/C
    disp('Please mind the step size. Numerical noice is coming... ')
        pause
end

%% pseudo-steady-state assumption
% advection
%   bulk temperature vectors
T_bulk_in_down=[T_in, T(1:end-1)];
T_bulk_in_up=[T(2:end),T_in];

%   mass flow rate vectors and output temperatures
m_dot_up_out=zeros(1,n_fv);
m_dot_down_out=zeros(1,n_fv);
if m_dot_top<0
    m_dot_up_out(1:i_cir)=abs(m_dot_top);
    outVars.T_top=T(1);
else
    m_dot_down_out(1:i_cir-1)=m_dot_top;
    outVars.T_top=T_in;
end
if m_dot_bot<0
    m_dot_up_out(i_cir+1:end)=abs(m_dot_bot);
    outVars.T_bot=T_in;
else
    m_dot_down_out(i_cir:end)=m_dot_bot;
    outVars.T_bot=T(end);
end
m_dot_up_in=[m_dot_up_out(2:end),m_dot_up_out(end)];
m_dot_down_in=[m_dot_down_out(1),m_dot_down_out(1:end-1)];

%   advection heat flow 
Q_dot_bulk=4187*(...
           -m_dot_down_out.*T...
           -m_dot_up_out.*T...
           +m_dot_down_in.*T_bulk_in_down...
           +m_dot_up_in.*T_bulk_in_up...
           );
Q_dot_bulk(i_cir)=Q_dot_bulk(i_cir)+4187*(m_dot_bot-m_dot_top)*T_in_cir;

% conduction
k=0.65; %W/mK thermal conductivity of water at 60°C
A=pi*(D/2)^2;%area of section between two finite volumes

T_cond_up=[T(1),T(1:end-1)];%no conduction at top
T_cond_low=[T(2:end),T(end)];%no conduction at bottum

Q_dot_con=k*A/(L/n_fv)*(T_cond_up-2*T+T_cond_low);

%convective losses to surroundings
Q_dot_loss=UA/n_fv*(T-T_env);%heat loss to surroundings

% heat transfer from spiral to tank
NTUi=UA_spi/(4187*m_dot_spi*n_fv);
T_term1=T_in_spi*exp(-[1:n_fv]*NTUi);
if m_dot_spi>0
    T_term2=(1-exp(-NTUi))*exp(-[1:n_fv]*NTUi).*cumsum(T.*exp([1:n_fv]*NTUi));
    T_spi_out=T_term1+T_term2;
    Q_dot_spi=-4187*m_dot_spi*diff([T_in_spi,T_spi_out]);
elseif m_dot_spi<0
    T_term2=(1-exp(-NTUi))*exp(-[1:n_fv]*NTUi).*cumsum(fliplr(T).*exp([1:n_fv]*NTUi));
    T_spi_out=T_term1+T_term2;
    Q_dot_spi=fliplr(-4187*m_dot_spi*diff([T_in_spi,T_spi_out]));
else
    Q_dot_spi=zeros(1,n_fv);% no transfer
    T_spi_out=nan(1,n_fv);
end

Q_dot=Q_dot_bulk... %energy flow due to bulk flow=advection
      + Q_dot_con... %conduction
      - Q_dot_loss... % losses to surroundings
      + Q_dot_spi; % heating by spiral
  
outVars.Q_dot_bulk=sum(Q_dot_bulk);
outVars.Q_dot_con=sum(Q_dot_con);
outVars.Q_dot_loss=sum(Q_dot_loss);
outVars.Q_dot_spi=sum(Q_dot_spi);
outVars.Q_dot=Q_dot;
outVars.T_spi_out=T_spi_out;

%% dynamics
outVars.T_fv=T+Q_dot*timestep/(4187*1000*L/n_fv*A);
% 
% if m_dot_top<0
%     keyboard
% end
end