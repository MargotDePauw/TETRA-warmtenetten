% testscript to test storage tank with variable step size 

timestep=60; %s
%% define params
n_fv=10;

param.D=0.5;%m
param.L=1;%m
param.UA=0.5*pi*param.D*param.L; %W/K
param.dT_max=1; %°C
param.n_TT=[1 n_fv];%[floor(1.1/10*n_fv) ceil(9.9/10*n_fv)]; %position of sensor
param.DT_hys=NaN;%°C
param.T_sp=[75 75];
param.T_env=20;
DTA=80-65;
DTB=60-10;
Q_dot_des=40000;
param.UA_spi=1.5*Q_dot_des/(DTA-DTB)/log(DTA/DTB); %W/K
m_dot_boi_nom=Q_dot_des/(4187*(80-60));

inPhys.T_sp=param.T_sp;

paramPipe.V=m_dot_boi_nom*timestep;
paramPipe.flag_Tknown=0;
paramPipe.UA=0;
inPhysPipe.T_env=20;

%% define physical inputs and timestep

timeTot=5*3600;  %total time in seconds
n_sim=timeTot/timestep; %number of simpulation points

T_in=10*ones(n_sim,1);
T_env=20*ones(n_sim,1);
m_dot=[0*ones(n_sim/2,1);-Q_dot_des/(4187*(60-10))*ones(n_sim/2,1)];

%% data storage and initialisation

T_out=zeros(n_sim,1);
T_spi_out=zeros(n_sim,n_fv);
T_spi_in=T_out;
T_fv=zeros(n_sim,n_fv);
Q_dot_loss=T_out;
Q_dot_spi=T_out;
e_loading=T_out;
m_dot_spi=T_out;

T_pipe_cont_mean=T_out;
T_pipe_out=T_out;

T_fv(1,:)=20*ones(n_fv,1);
T_spi_out(1)=20;
e_loading(1)=1;

inHisPipe.V=paramPipe.V;
inHisPipe.T=-1;


%% simulations


for i=2:n_sim
    % mass flow rates
    m_dot_boi=m_dot_boi_nom*e_loading(i-1);
    m_dot_spi(i)=m_dot_boi;
    inPhys.m_dot=m_dot(i);
    
    % return pipe (no supply pipe)
    inPhysPipe.T_in=T_spi_out(i-1);    
    inPhysPipe.m_dot=m_dot_boi;

    outPipe=pipeSingle_20180105(inPhysPipe,inHisPipe,...
                paramPipe,timestep);

    T_pipe_cont_mean(i-1)=outPipe.T_cont_mean;
    T_pipe_out(i)=outPipe.T_out;
    
    inHisPipe.T=outPipe.T_cont_new;
    inHisPipe.V=outPipe.V_cont_new;   
    
    % boiler
    if e_loading(i-1)
        Q_dot_boi=4187*m_dot_boi*(80-T_pipe_out(i));
        T_boi_out=T_pipe_out(i)+min(Q_dot_des,Q_dot_boi)/(4187*m_dot_boi);
    else
        Q_dot_boi=0;
        T_boi_out=T_pipe_out(i);
    end
    T_spi_in(i)=T_boi_out;
    
    %inputs
    inHis.T_fv=T_fv(i-1,:);
    inHis.e_loading=e_loading(i-1);
    inPhys.T_in=T_in(i);
    inPhys.T_in_spi=T_boi_out;
    inPhys.T_env=T_env(i);
    inPhys.m_dot_spi=m_dot_spi(i);
    
    %inPhys.T_sp=%defined at param list  
    %inPhysSto.T_env=20;%defined at param list    
    %call
    outVars=storageTankFV_varSS_loadCtr_Spiral_20180912(inPhys,inHis,param,timestep);
    %outputs
    T_out(i)=outVars.T_top;
    T_spi_out(i,:)=outVars.T_spi_out;
    T_fv(i,:)=outVars.T_fv;
    Q_dot_loss(i)=outVars.Q_dot_loss;
    Q_dot_spi(i)=outVars.Q_dot_spi;
    e_loading(i)=outVars.e_loading;
%     if i>=178
%         keyboard
%     end
end





%% plottings and analysis
% disp('total mass [kg], total time to flush tank [s] and current total time:...............  ')
mass_total=param.L*pi*(param.D/2)^2*1000; %in kg
time_flush=mass_total/1.2;
time_total=n_sim*timestep;

% plot normal
figure;
nsp=4; %number of suplots
t=timestep*(0:n_sim-1);
subplot(nsp,1,1)
    plot(t,T_fv,'-b')
    hold on
    plot(t,T_fv(:,1),'-r')
    plot(t,T_fv(:,end),'-r')
    plot(t,param.T_sp(1)*length(t),'-k')
    xlabel('time [s]'); ylabel('T [C]');
    ylim([0 100]);
    title('main content')
subplot(nsp,1,2)
    plot(t,T_spi_in,'-xr')
    hold on
    plot(t,T_spi_out,'-b')
    xlabel('time [s]'); ylabel('T [C]');
    ylim([0 100]);
    legend('in','out')
    title('spiral')
subplot(nsp,1,3)
    plot(t,m_dot,'-x')
    hold on
    plot(t,m_dot_spi,'-xr')
    xlabel('time [s]'); ylabel('m dot [kg/s]');% plot normal
    legend('main content','spiral')
    title('spiral')
subplot(nsp,1,4)
    plot(t,Q_dot_spi)
    xlabel('time [s]'); ylabel('Q dot [J/s]');% plot normal
    
