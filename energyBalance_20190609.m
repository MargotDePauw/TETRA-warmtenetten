% This function calculates the energy balances of a system of components or
% single components. It calculates all energy term for the given simulation
% results as input (number of points=n_sim), for a given time span.
%
% INPUTS
%   inHyd.c      [J/kgK] specific thermal capacity of fluid (normally
%                   water: 4187 and or air) ((vector) (n_simxn_in)
%   inHyd.m_dot_in[kg/s] mass flow rate of hydronic flow of inlet(s)(vector) (n_simxn_in)
%   inHyd.m_dot_out[kg/s] mass flow rate of hydronic flow of outlet(s)(vector) (n_simxn_out)
%   inHyd.T_in   [°C] temperature of ingoing water (vector) (n_simxn_in)
%   inHyd.T_out  [°C] temperature of outgoing water (vector) (n_simxn_out)
%
%   inHeatSource [W] a heat source applied to the system, for instance the
%                heat source in a boiler (n_simxn_heatSources)
%
%   inHeatLoss   [W] losses to surroundings (n_simxn_heatLosses)
%
%   inStor.C     [J/K] total thermal capacity of one ore more subsystems
%                with thermal inertia (1xn_capacities)
%   inStor.T     [°C] temperatures of upper capacities during simulation
%                (n_simxn_capacities).  
%   inRef.P      [W] reference energy flow: used to calculate relative
%                error. For heat production systems: take nominal heat
%                production; for storage system: take loading at nominal
%                rate; for building and rest of HVAC system: take nominal
%                heat emission by emitters
%
%   timestep     [s] (1x1)
%
% OUTPUTS
%
% outVars.tot ----- this gives values of the total period
%
%   E_hyd_in        [J] total ingoing hydronic/bulk flow energy (1x1)
%   E_hyd_out       [J] total outgoing hydronic/bulk flow energy (1x1)
%   E_hyd_net       [J] net hydronic energy flow (1x1)
%   E_heatSource    [J] total ingoing energy from heat source(s)(1x1)
%   E_heatLoss      [J] total outgoing energy from heat losses(1x1)
%   E_stor_net      [J] net total energy that is stored: end - begin (1x1)
%   E_InOut_net     [J] net total energy that goes in - out (1x1)
%   E_net           [J] total net energy (1x1):
%                   (IN-OUT) - STORAGE, should be (close to) zero
%   E_ref           [J] energy of reference
%   eRelHeatSource  [-] relative error: E_net/E_heatSource (1x1):
%                   how much energy 'disapears'
%   eRelHyd         [-] relative error: E_net/abs(E_hyd_in-E_hyd_out)(1x1):
%                   how much energy 'disapears'
%   eRelRef         [-] relative error: E_net/E_ref:
%                   how much energy 'disapears'
%
%
% outVars.instant ----- this gives values at moment i and of the period 
%                       between i-1 and i
%
%   P_hyd_in        [W] total hydronic heat flow in (n_sim,1)
%   P_hyd_out       [W] total hydronic heat flow out (n_sim,1)
%   P_heatSource    [W] total heat flow in from external sources (n_sim,1)
%   P_heatLoss      [W] total heat flow out to surroundings (n_sim,1)
%   E_hyd_in        [J] total ingoing hydronic/bulk flow energy (n_simx1)
%   E_hyd_out       [J] total outgoing hydronic/bulk flow energy (n_simx1)
%   E_hyd_net       [J] net hydronic energy flow (n_simx1)
%   E_heatSource    [J] total ingoing energy from heat source(s)(n_simx1)
%   E_heatLoss      [J] total outgoing energy from heat losses(n_simx1)
%   E_stor          [J] total energy that is stored (n_sim,1)
%   E_InOut_net     [J] net total energy that goes in - out (n_simx1)
%   E_net           [J] total net energy:
%                   (IN-OUT) - STORAGE, should be (close to) zero
%   E_ref           [J] energy of reference
%   eRelHeatSource  [-] relative error: E_net/E_heatSource:
%                   how much energy 'disapears'
%   eRelHyd         [-] relative error: E_net/abs(E_hyd_in-E_hyd_out):
%                   how much energy 'disapears'
%   eRelRef         [-] relative error: E_net/E_ref:
%                   how much energy 'disapears'
%
% outVars.summary ---- most important outputs, only these should be saved
%   eRelRef         see outVars.instant(n_simx1)
%   mean            mean of eRelRef (1x1)
%   median          median of eRelRef (1x1)
%   p10             10% percentile of eRelRef (1x1)
%   p25             25% percentile of eRelRef (1x1)
%   p75             75% percentile of eRelRef (1x1)
%   p90             90% percentile of eRelRef (1x1)
%   mean            mean of eRelRef (1x1)
%   min             min of eRelRef (1x1)
%   max             max of eRelRef (1x1)
%   eRelHeatSource  see outVars.tot (1x1)
%   eRelHyd         see outVars.tot (1x1)
%
% freek.vanriet@uantwerpen.be
% adjustments
% 2017 08 28    make 0xNaN equal to 0 for hydronic energy flows (bug fixed)
%               extra eRel added
%               E_stor_net erroneous calculation replaced
% 2017 08 30    outputs added and structure of output changed
% 2018 03 01    extra input: some reference energy flow to calculate
%               reference relative error
%               also a summary added which includes most important outputs
%               (partly same as other outputs)
%               also incremental data of previous version deleted
%               allow for multiple hydronic inlets (n_in) and outlets
%               (n_out)
% 2019 06 09    muliptle specific capacities for hydronics

function outVars=energyBalance_20190609(inHyd,inHeatSource,inHeatLoss,inStor,inRef,timestep)
%% hydronic energy
% make sure m_dot=0 results in P=0
inHyd.T_in(inHyd.m_dot_in==0)=0;
inHyd.T_out(inHyd.m_dot_out==0)=0;

% energy flows
P_hyd_in=sum(inHyd.c.*inHyd.m_dot_in.*inHyd.T_in,2);
P_hyd_out=sum(inHyd.c.*inHyd.m_dot_out.*inHyd.T_out,2);

% energy
%   total
E_hyd_in=timestep*sum(P_hyd_in);
E_hyd_out=timestep*sum(P_hyd_out);
E_hyd_net=E_hyd_in-E_hyd_out;

%   instataneous
E_hyd_in_instant=timestep*P_hyd_in;
E_hyd_out_instant=timestep*P_hyd_out;
E_hyd_net_instant=E_hyd_in_instant-E_hyd_out_instant;

%% heat source
% energy flow
P_heatSource=sum(inHeatSource,2); % according to second dimension (see dim)
% energy
%   total
E_heatSource=timestep*sum(P_heatSource);
%   instantaneous
E_heatSource_instant=timestep*P_heatSource;

%% heat losses

P_heatLoss=sum(inHeatLoss,2); % according to second dimension (see dim)

% energy
%   total
E_heatLoss=timestep*sum(P_heatLoss);
%   instantaneous
E_heatLoss_instant=timestep*P_heatLoss;

%% storage

% seperated enthalpy contents
[n_sim,~]=size(inStor.T);
C=repmat(inStor.C,n_sim,1);
E_stor=C.*inStor.T;

% total enthalpy contents
E_stor_tot=sum(E_stor,2); % total in the meaning of
                          % the sum of collumns
%   total (in time)
E_stor_net=E_stor_tot(end)-E_stor_tot(1);
%   instantaneous
E_stor_net_instant=diff(E_stor_tot);

%% reference
E_ref_instant=inRef.P*timestep;
E_ref_tot=sum(E_ref_instant);%was n_sim*E_ref_instant (MDP)

%% balances

% netto in -outgoing
%   total
E_InOut_net= E_hyd_net+E_heatSource-E_heatLoss;
%   instantaneous
E_InOut_net_instant=E_hyd_net_instant+E_heatSource_instant-E_heatLoss_instant;

% absolute balans
%   total
E_net=E_InOut_net-E_stor_net;
%   instantaneous
E_net_instant=E_InOut_net_instant(2:end)-E_stor_net_instant;

% relative errors
%   total
eRelHeatSource=abs(E_net/E_heatSource);
eRelHyd=abs(E_net/(E_hyd_in-E_hyd_out));
eRelRef=abs(E_net/E_ref_tot);
%   instantaneous
eRelHeatSource_instant=abs(E_net_instant./E_heatSource_instant(2:end));
eRelHyd_instant=abs(E_net_instant./(E_hyd_in_instant(2:end)-E_hyd_out_instant(2:end)));
eRelInOut=abs(E_net_instant./E_InOut_net_instant(2:end));
eRelRef_instant=abs(E_net_instant./E_ref_instant(2:end));%E_ref_instant(2:end) instead of E_ref_instant (MDP) and ./ instead of /

%% make summary
% outVars.summary ---- most important outputs, only these should be saved
outMean=mean(eRelRef_instant);
outMedian=median(eRelRef_instant);
outMin=min(eRelRef_instant);
outMax=max(eRelRef_instant);
% p10=prctile(eRelRef_instant,10);
% p25=prctile(eRelRef_instant,25);
% p75=prctile(eRelRef_instant,75);
% p90=prctile(eRelRef_instant,90);

%% assign outputs
% total
outVars.tot.E_hyd_in=E_hyd_in;
outVars.tot.E_hyd_out=E_hyd_out;
outVars.tot.E_hyd_net=E_hyd_net;
outVars.tot.E_heatSource=E_heatSource;
outVars.tot.E_heatLoss=E_heatLoss;
outVars.tot.E_stor_net=E_stor_net;
outVars.tot.E_InOut_net=E_InOut_net;
outVars.tot.E_net=E_net;
outVars.tot.E_ref=E_ref_tot;

outVars.tot.eRelHeatSource=eRelHeatSource;
outVars.tot.eRelHyd=eRelHyd;
outVars.tot.eRelRef=eRelRef;

% instantaneous

outVars.instant.E_hyd_in=E_hyd_in_instant;
outVars.instant.E_hyd_out=E_hyd_out_instant;
outVars.instant.E_hyd_net=E_hyd_net_instant;
outVars.instant.E_heatSource=E_heatSource_instant;
outVars.instant.E_heatLoss=E_heatLoss_instant;
outVars.instant.E_stor=E_stor;
outVars.instant.E_stor_tot=E_stor_tot;
outVars.instant.E_stor_net=E_stor_net_instant;
outVars.instant.E_InOut_net=E_InOut_net_instant;
outVars.instant.E_net=E_net_instant;
outVars.instant.E_ref=E_ref_instant;

outVars.instant.P_hyd_in=P_hyd_in;
outVars.instant.P_hyd_out=P_hyd_out;
outVars.instant.P_heatSource=P_heatSource;
outVars.instant.P_heatLoss=P_heatLoss;

outVars.instant.eRelHeatSource=eRelHeatSource_instant;
outVars.instant.eRelHyd=eRelHyd_instant;
outVars.instant.eRelInOut=eRelInOut;
outVars.instant.eRelRef=eRelRef_instant;


% summary
outVars.summary.mean=outMean;
outVars.summary.median=outMedian;
outVars.summary.min=outMin;
outVars.summary.max=outMax;
% outVars.summary.p10=p10;
% outVars.summary.p25=p25;
% outVars.summary.p75=p75;
% outVars.summary.p90=p90;
outVars.summary.eRelRef_instant=eRelRef_instant;
outVars.summary.eRelRef=eRelRef;
outVars.summary.eRelHeatSource=eRelHeatSource;
outVars.summary.eRelHyd=eRelHyd;

end



