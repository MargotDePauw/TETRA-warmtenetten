
n_sim=105120;%105120 for 5 minutes
timestep = 5*60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Het pad moet telkens worden aangepast!!
%pad ='C:\Users\u0035550\Documents\MATLAB 2017-2018\componentmodellen\sat unit\';% 'E:\Simulaties\Ecodroom klein individueel\standaard\'; 
%pad= 'C:\Users\u0035550\Thomas More\KCE - Team - Documenten\B1_Projecten-gebouwen\TET-2019-warmtenetten\model\model zonnecollector\matlab\';
pad='C:\Users\u0035550\Documents\TETRA warmtenetten\simulaties\model zonnecollector\code\';
save('pad','pad');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = [pad,'zonnedata_Freek.xlsx'];

weather = xlsread(file,'A3:J105123');

Te=weather(1:n_sim,2);
Sol_Zenith=weather(1:n_sim,4);
Sol_Azimuth=weather(1:n_sim,5);
I_beam=weather(1:n_sim,6)/3.6;%Watt/m²
I_dsky=weather(1:n_sim,7)/3.6;%Watt/m²
I_dgnd=weather(1:n_sim,8)/3.6;%Watt/m²
I_d=weather(1:n_sim,9)/3.6;%Watt/m²
AI=weather(1:n_sim,10);



%parameters for evacuated tubular collector CPC 1512 (fiche Sanutal_1)
% rend0=0.561, c1=0.45, c2=0.007
% flow rate for testing 0.02 kg/sm² (*2.28m²)
% C=78.478 kJ/Km²
param.A = 2.28;%m²
param.C=78.478*param.A*1000;  %J/°C
param.k0=0.561;
param.k1=0.45;
param.k2=0.007;
param.c_water= 4187; %J/kgK
param.m_dot_test=0.02*param.A;
param.beta_collect = 40; 
% manufacturer data IAMt and IAMl for transversal and longitudinal direction,column=theta_in
param.IAM.theta_t=[0 10 20 30 40 50 60 70 80 90];
param.IAM.IAM_t= [1 1.01 1.02 1.03 1.04 1.03 1.12 1.57 0.79 0.00];
param.IAM.theta_l=[0 10 20 30 40 50 60 70 80 90];
param.IAM.IAM_l= [1 0.98 0.96 0.96 0.96 0.88 0.78 0.68 0.34 0.00];
param.azimuth_collect=0;
param.ro_g= 0.2; %??

%input data
%inPhys.m_dot = 0.02*param.A;   %kg/s
inPhys.T_in = 20;



output = ones(n_sim,7);
output(1,:)= [20 0 0 0 0 0 0 ];

const = false;

for i=2:n_sim
    
inHis.T_out= output(i-1,1);
m_dot_prev = output(i-1,7);

%inPhys
if const == false
    inPhys.T_amb= Te(i);
    inPhys.I_b= I_beam(i);
    inPhys.I_dsky= I_dsky(i);
    inPhys.I_dgnd =I_dgnd(i);
    inPhys.azimuth_sol=Sol_Azimuth(i);
    inPhys.zenith_sol=Sol_Zenith(i);
    inPhys.theta=AI(i);
else
    inPhys.T_amb= (20+inHis.T_out)/2;
    inPhys.I_b= 2510/3.6;
    inPhys.I_dsky= 893/3.6;
    inPhys.I_dgnd =67/3.6;
    inPhys.azimuth_sol=35.3;
    inPhys.zenith_sol=39;
    inPhys.theta=22.3;    
end



if inHis.T_out > -9999 %25
    inPhys.m_dot = 0.02*param.A;   %kg/s
else
    if inHis.T_out<=22
        inPhys.m_dot = 0;   %kg/s
    else
        inPhys.m_dot=m_dot_prev;
    end
end
if i == 2796
    tijdstap=i;
end

outSolCol=solCol_20200528(inPhys,inHis,param,timestep);

output(i,:) = [outSolCol.T_out outSolCol.IAM  outSolCol.IAMb outSolCol.Q_dot outSolCol.Q_dot_loss outSolCol.corr_flow outSolCol.m_dot];
output_balance(i,:) = [outSolCol.T_out_mean outSolCol.m_dot outSolCol.Q_dot_rad outSolCol.Q_dot_loss outSolCol.Tcoll outSolCol.Q_dot outSolCol.Q_dot_rad outSolCol.DT_0];
end

%%

DT = output_balance(2:end,8);
ind1=(DT>5);
ind2=(DT<= 5);
l1=sum(ind1);
l2=sum(ind2);
ind=ind1;
%test energy balans
%vectorize constant input
n=n_sim-1; %l1; %n_sim
M_T_in = repmat(inPhys.T_in,n,1);
inHyd.c=param.c_water;
inHyd.m_dot_in = output_balance(2:end,2);%output_balance(ind,2); %kg/s  output_balance(:,2)
inHyd.m_dot_out = output_balance(2:end,2);%output_balance(ind,2); %kg/s  output_balance(:,2)
inHyd.T_in = M_T_in;
inHyd.T_out = output_balance(2:end,1);%output_balance(ind,1);%output_balance(:,1);
inHeatSource = output_balance(2:end,3);%output_balance(ind,3); %Watt  output_balance(:,3)
inHeatLoss = output_balance(2:end,4);%output_balance(ind,4); %Watt  output_balance(:,4)
inStor.T=output_balance(2:end,5);%output_balance(ind,5); %output_balance(:,5)
inStor.C= param.C; %J/K
%reference value
%inRef.P=output_balance(:,6); %Watt, netto energy gain, but not good as ref value 
Q_dot_rad = output_balance(2:end,7);%output_balance(ind,7); %solar gains to collector fluid output_balance(:,7)
Pref= max(Q_dot_rad); %max value of solar gains as reference value
inRef.P=repmat(Pref,n,1);
outBalanceCollector = energyBalance_20190609(inHyd,inHeatSource,inHeatLoss,inStor,inRef,timestep);










% figure
% ax(1)=subplot(2,1,1); hold on %secundaire kant
% plot(output(:,1),'-xr');
% plot(T22,'-xb');
% plot(T21,'xg');
% xlabel('time [sec]','FontSize',14)
% ylabel('TDHW[°C]','FontSize',14)
% h_legend=legend('Tcalc','Tmeasured','TC');
% set(h_legend,'FontSize',14);
% ax(2)=subplot(2,1,2); hold on
% plot(output(:,2),'-xr');
% plot(T13,'-xb');
% plot(T12,'-xg');
% xlabel('time [sec]','FontSize',14)
% ylabel('Tretour[°C]','FontSize',14)
% h_legend=legend('Tcalc','Tmeasured','Tsup');
% set(h_legend,'FontSize',14);
% linkaxes(ax,'x');
% zoom on;
% 
% 
% 
% figure
% ax(1)=subplot(4,1,1); hold on %secundaire kant
% plot(output(:,1),'-xr');
% plot(T22,'-xb');
% plot(T21,'xg');
% xlabel('time []','FontSize',14)
% ylabel('temperatures[°C]DHW','FontSize',14)
% h_legend=legend('Tcalc','TDHWm','TCm');
% set(h_legend,'FontSize',14);
% ax(2)=subplot(4,1,2); hold on
% plot(output(:,2),'-xr');
% plot(T13,'-xb');
% plot(T12,'-xg');
% xlabel('time []','FontSize',14)
% ylabel('temperatures[°C]retour','FontSize',14)
% h_legend=legend('Tcalc','Tretm','Tsupm');
% set(h_legend,'FontSize',14);
% ax(3)=subplot(4,1,3); hold on
% plot(output(:,3),'-xr');
% xlabel('time []','FontSize',14)
% ylabel('Tsatunit','FontSize',14)
% h_legend=legend('Tcalc');
% set(h_legend,'FontSize',14);
% ax(4)=subplot(4,1,4); hold on
% plot(FT12,'-xr');
% plot(FT21,'-xb');
% xlabel('time []','FontSize',14)
% ylabel('flow','FontSize',14)
% h_legend=legend('Fprim','Fsec');
% set(h_legend,'FontSize',14);
% linkaxes(ax,'x');
% zoom on;
% 
