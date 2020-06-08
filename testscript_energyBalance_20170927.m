% script to test energy Balances of simulation

% INPUTS defined acc. to simscript_boilerCHP_parallel_20170807
%% CHP
inHyd.c=4187;
inHyd.m_dot=storeOutCHP.m_dot_sec;
inHyd.T_in=storeOutPipeCHPr.T_out;
%inHyd.T_out=storeOutCHP.T_out;
inHyd.T_out=storeOutCHP.T_out_mean;
inHeatSource=storeOutCHP.Q_dot;
inHeatLoss=storeOutCHP.Q_dot_loss;
inStor.C=paramCHP.C;
inStor.T=storeOutCHP.T_out;
timestep=timestep;

outBalanceCHP=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

%% Boiler
inHyd.c=4187;
inHyd.m_dot=storeOutBoiler.m_dot;
inHyd.T_in=storeOutPipeZoneR.T_out;
%inHyd.T_out=storeOutBoiler.T_out;
inHyd.T_out=storeOutBoiler.T_out_mean;
inHeatSource=storeOutBoiler.PLR*paramBoiler.P_max;
inHeatLoss=storeOutBoiler.Q_dot_loss;
inStor.C=paramBoiler.C_boiler;
inStor.T=storeOutBoiler.T_out;
timestep=timestep;

outBalanceBoiler=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

%% Storage vessel
inHyd.c=4187;
inHyd.m_dot=storeOutStorage.m_dot;
inHyd.T_in=storeOutStorage.T_top; % ingoing temperature, if m_dot>0
inHyd.T_out=storeOutStorage.T_bot; % outgoing temperature, if m_dot>0
inHeatSource=zeros(n_sim,1);
inHeatLoss=storeOutStorage.Q_dot_loss;
C_storage=4187*1000*paramStorage.L/n_fv*pi*(paramStorage.D/2)^2*ones(1,n_fv);
inStor.C=C_storage;
inStor.T=storeOutStorage.T_fv;
timestep=timestep;

outBalanceStorage=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

%% Zone

% bounderies emitter
inHyd.c=4187;
inHyd.m_dot=storeOutDistrFlow.m_dot;
inHyd.T_in=storeOutPipeZoneS.T_out; 
%inHyd.T_out=storeOutZone.T_out; 
inHyd.T_out=storeOutZone.T_out_mean; 
inHeatSource=zeros(n_sim,1);
inHeatLoss=storeOutZone.P_em;
inStor.C=[paramZone.C_em];
inStor.T=[storeOutZone.T_out];
timestep=timestep;

outBalanceZone1=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

% bounderies air
inHyd.c=4187;
inHyd.m_dot=storeOutDistrFlow.m_dot;
inHyd.T_in=storeOutPipeZoneS.T_out; 
%inHyd.T_out=storeOutZone.T_out; 
inHyd.T_out=storeOutZone.T_out_mean; 
inHeatSource=zeros(n_sim,1);
inHeatLoss=storeOutZone.P_em_wall+storeOutZone.P_zone_wall;
inStor.C=[paramZone.C_em,paramZone.C_zone];
inStor.T=[storeOutZone.T_out,storeOutZone.T_zone];
timestep=timestep;

outBalanceZone2=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

% bounderies total zone
inHyd.c=4187;
inHyd.m_dot=storeOutDistrFlow.m_dot;
inHyd.T_in=storeOutPipeZoneS.T_out; 
%inHyd.T_out=storeOutZone.T_out; 
inHyd.T_out=storeOutZone.T_out_mean;  
inHeatSource=zeros(n_sim,1);
inHeatLoss=storeOutZone.P_wall_ext;
inStor.C=[paramZone.C_em,paramZone.C_zone,paramZone.C_wall];
inStor.T=[storeOutZone.T_out,storeOutZone.T_zone,storeOutZone.T_wall];
timestep=timestep;

outBalanceZone3=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

%% Pipes

% Return, to CHP
inHyd.c=4187;
inHyd.m_dot=storeOutCHP.m_dot_sec;
inHyd.T_in=storeOutTpiece1.T_out_return; 
inHyd.T_out=storeOutPipeCHPr.T_out;
inHeatSource=zeros(n_sim,1);
inHeatLoss=zeros(n_sim,1);
inStor1.C=inHyd.c*paramPipeCHPr.V;

maxLength=max(cellfun(@(x)numel(x),storeOutPipeCHPr.T_cont));
T_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeCHPr.T_cont,'UniformOutput',false));
V_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeCHPr.V_cont,'UniformOutput',false));
inStor1.T=sum(T_cont_matrix.*V_cont_matrix,2)/paramPipeCHPr.V;

timestep=timestep;

outBalancePipeCHPr=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor1,timestep);

% Supply, from CHP
inHyd.c=4187;
inHyd.m_dot=storeOutCHP.m_dot_sec;
%inHyd.T_in=storeOutCHP.T_out; 
inHyd.T_in=storeOutCHP.T_out_mean; 
inHyd.T_out=storeOutPipeCHPs.T_out; 
inHeatSource=zeros(n_sim,1);
inHeatLoss=zeros(n_sim,1);
inStor2.C=inHyd.c*paramPipeCHPs.V;

maxLength=max(cellfun(@(x)numel(x),storeOutPipeCHPs.T_cont));
T_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeCHPs.T_cont,'UniformOutput',false));
V_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeCHPs.V_cont,'UniformOutput',false));
inStor2.T=sum(T_cont_matrix.*V_cont_matrix,2)/paramPipeCHPs.V;

timestep=timestep;

outBalancePipeCHPs=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor2,timestep);

% Return, from Zone
inHyd.c=4187;
inHyd.m_dot=storeOutDistrFlow.m_dot;
% inHyd.T_in=storeOutZone.T_out; 
inHyd.T_in=storeOutZone.T_out_mean; 
inHyd.T_out=storeOutPipeZoneR.T_out; 
inHeatSource=zeros(n_sim,1);
inHeatLoss=zeros(n_sim,1);
inStor3.C=inHyd.c*paramPipeZoneR.V;

maxLength=max(cellfun(@(x)numel(x),storeOutPipeZoneR.T_cont));
T_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeZoneR.T_cont,'UniformOutput',false));
V_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeZoneR.V_cont,'UniformOutput',false));
inStor3.T=sum(T_cont_matrix.*V_cont_matrix,2)/paramPipeZoneR.V;

timestep=timestep;

outBalancePipeZoneR=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor3,timestep);

% Supply, to zone
inHyd.c=4187;
inHyd.m_dot=storeOutDistrFlow.m_dot;
inHyd.T_in=storeOutTpiece3.T_out; 
inHyd.T_out=storeOutPipeZoneS.T_out; 
inHeatSource=zeros(n_sim,1);
inHeatLoss=zeros(n_sim,1);
inStor4.C=inHyd.c*paramPipeZoneS.V;

maxLength=max(cellfun(@(x)numel(x),storeOutPipeZoneS.T_cont));
T_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeZoneS.T_cont,'UniformOutput',false));
V_cont_matrix=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),storeOutPipeZoneS.V_cont,'UniformOutput',false));
inStor4.T=sum(T_cont_matrix.*V_cont_matrix,2)/paramPipeZoneS.V;

timestep=timestep;

outBalancePipeZoneS=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor4,timestep);

%% total system
inHyd.c=4187;
inHyd.m_dot=zeros(n_sim,1);
inHyd.T_in=zeros(n_sim,1);
inHyd.T_out=zeros(n_sim,1);
inHeatSource=[storeOutCHP.Q_dot, storeOutBoiler.PLR*paramBoiler.P_max];
inHeatLoss=[storeOutCHP.Q_dot_loss, storeOutBoiler.Q_dot_loss,storeOutStorage.Q_dot_loss, storeOutZone.P_wall_ext];
inStor.C=[paramCHP.C,...
          paramBoiler.C_boiler,...
          C_storage,...
          paramZone.C_em, paramZone.C_zone, paramZone.C_wall,...
          inStor1.C,inStor2.C,inStor3.C,inStor4.C];
inStor.T=[storeOutCHP.T_out,...
          storeOutBoiler.T_out, ...
          storeOutStorage.T_fv, ...
          storeOutZone.T_out, storeOutZone.T_zone, storeOutZone.T_wall, ...
          inStor1.T,inStor2.T,inStor3.T,inStor4.T];  

timestep=timestep;

outBalanceTot=energyBalance_20170830(inHyd,inHeatSource,inHeatLoss,inStor,timestep);

%% plotting 1: absolute errors
% % CHP
% 
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: CHP')
%     plot(t(2:end),outBalanceCHP.instant.E_net,'-x');   
%     ylabel('error [J]','FontSize',14)
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalanceCHP.instant.E_stor_net,'-x');  
%     plot(t,outBalanceCHP.instant.E_hyd_net,'-x');
%     plot(t,outBalanceCHP.instant.E_heatSource,'-x');
%     plot(t,outBalanceCHP.instant.E_heatLoss,'-x');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
%     plot(t,outBalanceCHP.instant.E_InOut_net,'-xb');   
%     plot(t(2:end),outBalanceCHP.instant.E_stor_net,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('IN - OUT','STORAGE')
%     
% % Boiler
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: Boiler')
%     plot(t(2:end),outBalanceBoiler.instant.E_net,'-x');   
%     ylabel('error [J]','FontSize',14)
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalanceBoiler.instant.E_stor_net,'-x');  
%     plot(t,outBalanceBoiler.instant.E_hyd_net,'-x');
%     plot(t,outBalanceBoiler.instant.E_heatSource,'-x');
%     plot(t,outBalanceBoiler.instant.E_heatLoss,'-x');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
%     plot(t,outBalanceBoiler.instant.E_InOut_net,'-xb');   
%     plot(t(2:end),outBalanceBoiler.instant.E_stor_net,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('IN - OUT','STORAGE')
%     
% % storage
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: vessel')
%     plot(t(2:end),outBalanceStorage.instant.E_net,'-x');   
%     ylabel('error [J]','FontSize',14)
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalanceStorage.instant.E_stor_net,'-x');  
%     plot(t,outBalanceStorage.instant.E_hyd_net,'-x');
%     plot(t,outBalanceStorage.instant.E_heatSource,'-x');
%     plot(t,outBalanceStorage.instant.E_heatLoss,'-x');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
%     plot(t,outBalanceStorage.instant.E_InOut_net,'-xb');   
%     plot(t(2:end),outBalanceStorage.instant.E_stor_net,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('IN - OUT','STORAGE')
% 
% 
% % zone (total zone sustem)
% 
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: zone + emitter')
%     plot(t(2:end),outBalanceZone3.instant.E_net,'-x');   
%     ylabel('error [J]','FontSize',14)
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalanceZone3.instant.E_stor_net,'-x');  
%     plot(t,outBalanceZone3.instant.E_hyd_net,'-x');
%     plot(t,outBalanceZone3.instant.E_heatSource,'-x');
%     plot(t,outBalanceZone3.instant.E_heatLoss,'-x');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
%     plot(t,outBalanceZone3.instant.E_InOut_net,'-xb');   
%     plot(t(2:end),outBalanceZone3.instant.E_stor_net,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('IN - OUT','STORAGE')
% 
% 
% 
% % pipes related to zone
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: pipes related to zone')
%     plot(t(2:end),outBalancePipeZoneR.instant.E_net,'-xb');   
%     hold on
%     plot(t(2:end),outBalancePipeZoneS.instant.E_net,'-xr');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('error [J]','FontSize',14)
%     legend('return pipe','supply pipe')
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalancePipeZoneR.instant.E_stor_net,'-xb');  
%     hold on
%     plot(t,outBalancePipeZoneR.instant.E_hyd_net,'-ob');
%     plot(t(2:end),outBalancePipeZoneS.instant.E_stor_net,'-xr');  
%     plot(t,outBalancePipeZoneS.instant.E_hyd_net,'-or');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('return pipe STORAGE','return pipe HYDRON','supply pipe STORAGE','supply pipe HYDRON')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
% title('mean temperature content')
%     plot(t,inStor3.T,'-xb');   
%     hold on
%     plot(t,inStor4.T,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('T [C]','FontSize',14)
%     legend('return pipe','supply pipe')    
%     
% % pipes related to CHP
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: pipes related to CHP')
%     plot(t(2:end),outBalancePipeCHPr.instant.E_net,'-xb');   
%     hold on
%     plot(t(2:end),outBalancePipeCHPs.instant.E_net,'-xr');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('error [J]','FontSize',14)
%     legend('return pipe','supply pipe')
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalancePipeCHPr.instant.E_stor_net,'-xb');  
%     hold on
%     plot(t,outBalancePipeCHPr.instant.E_hyd_net,'-ob');
%     plot(t(2:end),outBalancePipeCHPs.instant.E_stor_net,'-xr');  
%     plot(t,outBalancePipeCHPs.instant.E_hyd_net,'-or');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('return pipe STORAGE','return pipe HYDRON','supply pipe STORAGE','supply pipe HYDRON')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
% title('mean temperature content')
%     plot(t,inStor1.T,'-xb');   
%     hold on
%     plot(t,inStor2.T,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('T [C]','FontSize',14)
%     legend('return pipe','supply pipe')  
% 
% 
% % total: HVAC + building
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: total system = HVAC + building')
%     plot(t(2:end),outBalanceTot.instant.E_net,'-x');   
%     ylabel('error [J]','FontSize',14)
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalanceTot.instant.E_stor_net,'-x');  
%     plot(t,outBalanceTot.instant.E_hyd_net,'-x');
%     plot(t,outBalanceTot.instant.E_heatSource,'-x');
%     plot(t,outBalanceTot.instant.E_heatLoss,'-x');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
%     plot(t,outBalanceTot.instant.E_InOut_net,'-xb');   
%     plot(t(2:end),outBalanceTot.instant.E_stor_net,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('IN - OUT','STORAGE')
% 
% 
% 
%  linkaxes(ax,'x');    
% %   

%% plotting 2: relative errors
n_ax=1;
t=timestep*[0:n_sim-1];
% CHP

figure;

ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
title('balance errors: CHP')
    plot(t(2:end),outBalanceCHP.instant.eRelInOut,'-x');   
    ylabel('rel error to in-out [-]','FontSize',14)
ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
    plot(t(2:end),outBalanceCHP.instant.E_stor_net,'-x');  
    plot(t,outBalanceCHP.instant.E_hyd_net,'-x');
    plot(t,outBalanceCHP.instant.E_heatSource,'-x');
    plot(t,outBalanceCHP.instant.E_heatLoss,'-x');
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
    plot(t,outBalanceCHP.instant.E_InOut_net,'-xb');   
    plot(t(2:end),outBalanceCHP.instant.E_stor_net,'-or');   
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('IN - OUT','STORAGE')
    
% Boiler
figure;

ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
title('balance errors: Boiler')
    plot(t(2:end),outBalanceBoiler.instant.eRelInOut,'-x');   
    ylabel('rel error to in-out [-]','FontSize',14)
ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
    plot(t(2:end),outBalanceBoiler.instant.E_stor_net,'-x');  
    plot(t,outBalanceBoiler.instant.E_hyd_net,'-x');
    plot(t,outBalanceBoiler.instant.E_heatSource,'-x');
    plot(t,outBalanceBoiler.instant.E_heatLoss,'-x');
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
    plot(t,outBalanceBoiler.instant.E_InOut_net,'-xb');   
    plot(t(2:end),outBalanceBoiler.instant.E_stor_net,'-or');   
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('IN - OUT','STORAGE')
    
% storage
figure;

ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
title('balance errors: vessel')
    plot(t(2:end),outBalanceStorage.instant.eRelInOut,'-x');   
    ylabel('rel error to in-out [-]','FontSize',14)
ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
    plot(t(2:end),outBalanceStorage.instant.E_stor_net,'-x');  
    plot(t,outBalanceStorage.instant.E_hyd_net,'-x');
    plot(t,outBalanceStorage.instant.E_heatSource,'-x');
    plot(t,outBalanceStorage.instant.E_heatLoss,'-x');
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
    plot(t,outBalanceStorage.instant.E_InOut_net,'-xb');   
    plot(t(2:end),outBalanceStorage.instant.E_stor_net,'-or');   
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('IN - OUT','STORAGE')


% zone (total zone sustem)

figure;

ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
title('balance errors: zone + emitter')
    plot(t(2:end),outBalanceZone3.instant.eRelInOut,'-x');   
    ylabel('rel error to in-out [-]','FontSize',14)
ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
    plot(t(2:end),outBalanceZone3.instant.E_stor_net,'-x');  
    plot(t,outBalanceZone3.instant.E_hyd_net,'-x');
    plot(t,outBalanceZone3.instant.E_heatSource,'-x');
    plot(t,outBalanceZone3.instant.E_heatLoss,'-x');
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
    plot(t,outBalanceZone3.instant.E_InOut_net,'-xb');   
    plot(t(2:end),outBalanceZone3.instant.E_stor_net,'-or');   
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('IN - OUT','STORAGE')



% % pipes related to zone
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: pipes related to zone')
%     plot(t(2:end),outBalancePipeZoneR.instant.eRelInOut,'-xb');   
%     hold on
%     plot(t(2:end),outBalancePipeZoneS.instant.eRelInOut,'-xr');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('rel error to in-out [-]','FontSize',14)
%     legend('return pipe','supply pipe')
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalancePipeZoneR.instant.E_stor_net,'-xb');  
%     hold on
%     plot(t,outBalancePipeZoneR.instant.E_hyd_net,'-ob');
%     plot(t(2:end),outBalancePipeZoneS.instant.E_stor_net,'-xr');  
%     plot(t,outBalancePipeZoneS.instant.E_hyd_net,'-or');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('return pipe STORAGE','return pipe HYDRON','supply pipe STORAGE','supply pipe HYDRON')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
% title('mean temperature content')
%     plot(t,inStor3.T,'-xb');   
%     hold on
%     plot(t,inStor4.T,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('T [C]','FontSize',14)
%     legend('return pipe','supply pipe')    
%     
% % pipes related to CHP
% figure;
% 
% ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
% title('balance errors: pipes related to CHP')
%     plot(t(2:end),outBalancePipeCHPr.instant.eRelInOut,'-xb');   
%     hold on
%     plot(t(2:end),outBalancePipeCHPs.instant.eRelInOut,'-xr');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('rel error to in-out [-]','FontSize',14)
%     legend('return pipe','supply pipe')
% ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
%     plot(t(2:end),outBalancePipeCHPr.instant.E_stor_net,'-xb');  
%     hold on
%     plot(t,outBalancePipeCHPr.instant.E_hyd_net,'-ob');
%     plot(t(2:end),outBalancePipeCHPs.instant.E_stor_net,'-xr');  
%     plot(t,outBalancePipeCHPs.instant.E_hyd_net,'-or');
%     xlabel('time [s]','FontSize',14)
%     ylabel('Energy [J]','FontSize',14)
%     legend('return pipe STORAGE','return pipe HYDRON','supply pipe STORAGE','supply pipe HYDRON')    
% ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
% title('mean temperature content')
%     plot(t,inStor1.T,'-xb');   
%     hold on
%     plot(t,inStor2.T,'-or');   
%     xlabel('time [s]','FontSize',14)
%     ylabel('T [C]','FontSize',14)
%     legend('return pipe','supply pipe')  


% total: HVAC + building
figure;

ax(n_ax)=subplot(3,1,1); n_ax=n_ax+1;hold on
title('balance errors: total system = HVAC + building')
    plot(t(2:end),outBalanceTot.instant.eRelInOut,'-x');   
    ylabel('rel error to in-out [-]','FontSize',14)
ax(n_ax)=subplot(3,1,2); n_ax=n_ax+1;hold on
    plot(t(2:end),outBalanceTot.instant.E_stor_net,'-x');  
    plot(t,outBalanceTot.instant.E_hyd_net,'-x');
    plot(t,outBalanceTot.instant.E_heatSource,'-x');
    plot(t,outBalanceTot.instant.E_heatLoss,'-x');
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('STORAGE','HYDRON','HEAT SOURCE','HEAT LOSS')    
ax(n_ax)=subplot(3,1,3); n_ax=n_ax+1;hold on
    plot(t,outBalanceTot.instant.E_InOut_net,'-xb');   
    plot(t(2:end),outBalanceTot.instant.E_stor_net,'-or');   
    xlabel('time [s]','FontSize',14)
    ylabel('Energy [J]','FontSize',14)
    legend('IN - OUT','STORAGE')



 linkaxes(ax,'x');    
%   