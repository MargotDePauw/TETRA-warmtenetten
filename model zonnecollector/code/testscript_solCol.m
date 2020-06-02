
n_sim=50000;
timestep = 5*60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Het pad moet telkens worden aangepast!!
%pad ='C:\Users\u0035550\Documents\MATLAB 2017-2018\componentmodellen\sat unit\';% 'E:\Simulaties\Ecodroom klein individueel\standaard\'; 
pad= 'C:\Users\u0035550\Thomas More\KCE - Team - Documenten\B1_Projecten-gebouwen\TET-2019-warmtenetten\model\model zonnecollector\matlab\';
save('pad','pad');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = [pad,'zonnedata_Freek.xlsx'];

weather = xlsread(file,'A3:J105123');

Te=weather(1:n_sim,2);
Sol_Zenith=weather(1:n_sim,4);
Sol_Azimuth=weather(1:n_sim,5);
I_beam=weather(1:n_sim,6)/3.6;%Watt
I_dsky=weather(1:n_sim,7)/3.6;%Watt
I_dgnd=weather(1:n_sim,8)/3.6;%Watt
I_d=weather(1:n_sim,9)/3.6;%Watt
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
param.beta_collect = 40; 
% manufacturer data IAMt and IAMl for transversal and longitudinal direction,column=theta_in
param.IAM.theta_t=[0 10 20 30 40 50 60 70 80 90];
param.IAM.IAM_t= [1 1.01 1.02 1.03 1.04 1.03 1.12 1.57 0.79 0.00];
param.IAM.theta_l=[0 10 20 30 40 50 60 70 80 90];
param.IAM.IAM_l= [1 0.98 0.96 0.96 0.96 0.88 0.78 0.68 0.34 0.00];
param.azimuth_collect=0;
param.ro_g= 0.2; %??

%input data
inPhys.m_dot = 0.02*param.A;   %kg/s
inPhys.T_in = 20;



output = ones(n_sim,6);
output(1,:)= [50 0 0 0 0 0];

for i=2:n_sim
%inPhys

inPhys.T_amb= Te(i);
inPhys.I_b= I_beam(i);
inPhys.I_dsky= I_dsky(i);
inPhys.I_dgnd =I_dgnd(i);
inPhys.azimuth_sol=Sol_Azimuth(i);
inPhys.zenith_sol=Sol_Zenith(i);
inPhys.theta=AI(i);

inHis.T_out= output(i-1,1);
inHis.t_run= output(i-1,2);


if i == 31803
    tijdstap=i;
end

outSolCol=solCol_20200528(inPhys,inHis,param,timestep);

output(i,:) = [outSolCol.T_out outSolCol.t_run  outSolCol.IAM  outSolCol.IAMb outSolCol.Q_dot_con outSolCol.Q_dot_loss];
end












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
