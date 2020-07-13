%% start
n_sim = 2000;
timestep = 5;
output = ones(n_sim,2);
Tboiler_init=10;
n_volumes=50;
outStorageInHis=Tboiler_init*ones(1,n_volumes);
output(1,:)= [1 20];

boiler=2; % 1=tesis BB, 2=collindi90L
% constant input/parameters
switch boiler
    case 1
       % param
       % storage tank
        param.D = 0.544;
        param.L=1.325;
        param.UA=0.74;%W/K
        param.dT_max=10;%????????????
        param.n_TT = 1;
        param.T_sp=1000;
        param.DT_hys = 1;
        %coil
        param.n_HE=[30 45];% position of HE (entering and leaving volumina)
        param.R_HE=(1/341)*16;%for entire HE/for 1 loop  *16

        inPhys.T_in= 10; %supply tank
        inPhys.T_env = 20;
        inPhys.m_dot=0;%water through tank
        inPhys.T_in_HE = 57; % temp of fluid entering heat exchanger1
        inPhys.m_dot_HE = 0.018; %flow heat exchanger1
    case 2
       % param
       % storage tank
        param.D = 0.4;
        param.L=0.72;
        param.UA=0.69;%W/K
        param.dT_max=10;%????????????
        param.n_TT = 1;
        param.T_sp=1000;
        param.DT_hys = 1;
        %coil
        param.n_HE=[25 50];% position of HE (entering and leaving volumina)
        param.R_HE=(1/483)*(max(param.n_HE)-min(param.n_HE)+1);%for entire HE/for 1 loop  /16

        inPhys.T_in= 10; %supply tank
        inPhys.T_env = 20;
        inPhys.m_dot=0;%water through tank
        inPhys.T_in_HE = 60; % temp of fluid entering heat exchanger1
        inPhys.m_dot_HE = 0.125; %flow heat exchanger1
end




for i=2:n_sim

    inHis.T_fv = outStorageInHis; %states of the finite volume of the storage tank
    inHis.e_loading=output(i-1,1);

    outStorage=storage_1coil(inPhys,inHis,param,timestep);
    output(i,:) = [outStorage.e_loading outStorage.T_ret];
    outStorageInHis=outStorage.T_fv;

    if rem(i,120)==0
        r=i/120;
        tempprofile (r,:)=outStorage.T_fv;
    end

end

