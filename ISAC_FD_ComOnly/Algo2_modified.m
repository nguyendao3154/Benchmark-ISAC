%*********************%
%  SYSTEM PARAMETERS  %
%*********************%
%% Some constants
rng('default');                     % Set random number generator for reproducibility
fc = 39e9;                          % Carrier frequency (Hz)
c = 3e8;                            % Speed of light
Nt =  32;                              % Number of array elements
Nr = 32;
%%
sUCA = antenna_geometry(Nt, fc);
steervector = phased.SteeringVector('IncludeElementResponse', true, "SensorArray",sUCA);
%% Some more fixed parameters
OFDM.N = 100;              % Number of subcarriers
OFDM.delta_f = 0.1e9;       % Carrier spacing
OFDM.carrier_freq = fc:OFDM.delta_f:fc+OFDM.N*OFDM.delta_f;
OFDM.lambda = c./OFDM.carrier_freq;
% pattern(sUCA,fc,0,[0:90],'Type','powerdb','PropagationSpeed',c)
gain = phased.ArrayGain('SensorArray',sUCA);
g = gain(fc,[zeros(1,91);0:90])
plot(1:91,db2mag(g));
%% Downlink
DL.NumUsers = 2;                              % Number of communication users
DL.NumSymbols = 1;                             % Number of communication symbols
DL.Constellation = 4;                              % QPSK
DL.data = randi([0 DL.Constellation-1], DL.NumUsers, DL.NumSymbols);        % Binary data
DL.S = pskmod(DL.data, DL.Constellation, pi/DL.Constellation);          % QPSK symbols
% DL.UsersPos = [rand(1, DL.NumUsers)*15; zeros(1, DL.NumUsers)*24; rand(1, DL.NumUsers)*10];
% DL.azimuth = atand(DL.UsersPos(1,:)./DL.UsersPos(2,:));
% DL.elevation = atand(DL.UsersPos(3,:)./sqrt(DL.UsersPos(1,:).^2 + DL.UsersPos(2,:).^2));
DL.range = 200*ones(1, DL.NumUsers);
% DL.azimuth = [0];
% DL.elevation = [70];
DL.azimuth = [0, 120];
DL.elevation = [80, 80];
DL.UsersPos = DL.range.*[cosd(DL.elevation).*sind(DL.azimuth); cosd(DL.elevation).*cosd(DL.azimuth); sind(DL.elevation)];
DL.sv = steervector(fc, [DL.azimuth; DL.elevation]);
DL.sv = DL.sv/sqrt(Nt);
% Channel
    DL.Pathloss = (OFDM.lambda(1)/(4*pi)/DL.range(1))^2;
    DL.channel = ones(1,DL.NumUsers)*sqrt(1e-12).*DL.sv/sqrt(1e-13);
    %% Uplink
UL.NumUsers = 2;                              % Number of communication users
UL.range = 200*ones(1, DL.NumUsers);
% UL.azimuth = [0];
% UL.elevation = [50];
UL.azimuth = [90, 180];
UL.elevation = [60, 60];
UL.UsersPos = DL.range.*[cosd(UL.elevation).*sind(UL.azimuth); cosd(UL.elevation).*cosd(UL.azimuth); sind(UL.elevation)];
UL.sv = steervector(fc, [UL.azimuth; UL.elevation]);
UL.sv = UL.sv./sqrt(Nr);

% Channel
    UL.Pathloss = OFDM.lambda(1)/(4*pi)/UL.range(1);
    UL.channel = ones(1,UL.NumUsers)*sqrt(1e-12).*UL.sv/sqrt(1e-13);
%% Calculate steering vector
ang_elv = linspace(0, 90, 181);       % Grid of Elevation angles
ang_az = linspace(-180, 180, 181);       % Grid of Elevation angles

sv_azimuth.DL1 = steervector(fc, [ang_az; DL.elevation(1).*ones(size(ang_az))]);
sv_elevation.DL1 = steervector(fc, [DL.azimuth(1).*ones(size(ang_elv)); ang_elv]);

sv_azimuth.DL2 = steervector(fc, [ang_az; DL.elevation(2).*ones(size(ang_az))]);
sv_elevation.DL2 = steervector(fc, [DL.azimuth(2).*ones(size(ang_elv)); ang_elv]);

sv_azimuth.UL1 = steervector(fc, [ang_az; UL.elevation(1).*ones(size(ang_az))]);
sv_elevation.UL1 = steervector(fc, [UL.azimuth(1).*ones(size(ang_elv)); ang_elv]);

sv_azimuth.UL2 = steervector(fc, [ang_az; UL.elevation(2).*ones(size(ang_az))]);
sv_elevation.UL2 = steervector(fc, [UL.azimuth(2).*ones(size(ang_elv)); ang_elv]);
% Sefl Interference (Rician)
K_rician_dB = 1;
K_rician = 10.^(K_rician_dB/10);
H_ray = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2); %CN(0,1)
H_los = ones(Nr,Nt);

%/Self interference cancellation level:
deta_dB = -90;
deta =   10.^(deta_dB/10);

Hli =(sqrt(K_rician*deta/(K_rician+1))*H_los+sqrt(deta/(K_rician+1))*H_ray);

%% Data scaling
NoisePower_dl = 1;
NoisePower_ul = 1;
scale_dl = NoisePower_dl;
DL.channel = DL.channel/sqrt(scale_dl); % scale DL channels with scaling factor of scale_dl
EffNoisePower_dl = NoisePower_dl/scale_dl; % the effective noise power after DL channels scaling

scale_ul = NoisePower_ul; 
UL.channel = UL.channel/sqrt(scale_ul); % scale UL channels with scaling factor of scale_ul 
Hli = Hli/sqrt(scale_ul); % scale SI channel with scaling factor of scale_ul
EffNoisePower_ul = NoisePower_ul/scale_ul; % effective noise power after UL channels scaling

%% generate initial points
Pu_W = 3.*ones(UL.NumUsers,1);
Pd_W = 64;
q_ul_init = Pu_W;
g_ul_dl = zeros(UL.NumUsers, DL.NumUsers);

Q_dl_init_tmp = zeros(Nt,Nt,DL.NumUsers);     % beamformer
for iUser_dl=1:DL.NumUsers
    a = rand(Nt,Nt)+1i*rand(Nt,Nt);
    Q_dl_init_tmp(:,:,iUser_dl) = a*a';
end
P_tmp = real(trace(sum(Q_dl_init_tmp,3)));
Q_dl_init = Q_dl_init_tmp/(P_tmp/(Pd_W));

mypsi = zeros(DL.NumUsers,1);
SINR_dl = zeros(DL.NumUsers,1);
t_dl_init = zeros(DL.NumUsers,1);
for iUser_dl=1:DL.NumUsers
    interf_dl = EffNoisePower_dl;
    interf_dl = interf_dl + real(DL.channel(:,iUser_dl)'*sum(Q_dl_init(:,:,1:DL.NumUsers~=iUser_dl),3)*DL.channel(:,iUser_dl));
    interf_dl = interf_dl + sum(q_ul_init.*abs(g_ul_dl(:,iUser_dl)).^2);
    SINR_dl(iUser_dl) = real((DL.channel(:,iUser_dl)'*Q_dl_init(:,:,iUser_dl)*DL.channel(:,iUser_dl)))/interf_dl;
    t_dl_init(iUser_dl) = 1+SINR_dl(iUser_dl);
    mypsi(iUser_dl) = (1+SINR_dl(iUser_dl))/interf_dl;
end
self_interference = EffNoisePower_ul*eye(Nr,Nr)+Hli*sum(Q_dl_init,3)*Hli';
X_init = zeros(Nr,Nr,UL.NumUsers);
x_ul_init = zeros(UL.NumUsers,1);
for iUser_ul=1:UL.NumUsers

    %X_init(:,:,iUser_ul)= interf_ul;
    X_init(:,:,iUser_ul) = self_interference+UL.channel(:,iUser_ul+1:UL.NumUsers)...
        *diag(q_ul_init(iUser_ul+1:UL.NumUsers))*UL.channel(:,iUser_ul+1:UL.NumUsers)';
    x_ul_init(iUser_ul)= sqrt(q_ul_init(iUser_ul));
end

ops=sdpsettings('solver','mosek','verbose',0); % solver options
% Define optimization variables
t_dl = sdpvar(DL.NumUsers,1);
t_ul = sdpvar(UL.NumUsers,1);
q_ul = sdpvar(UL.NumUsers,1); % uplink power 
Q_dl = sdpvar(Nt,Nt,DL.NumUsers,'hermitian','complex'); % downlink covariances
mybeta = sdpvar(DL.NumUsers,1) ;
x_ul = sdpvar(UL.NumUsers,1);
X = zeros(Nr,Nr,UL.NumUsers);

% define (fixed) constraints:
F=[];
for iUser_dl=1:DL.NumUsers
    F=[F,Q_dl(:,:,iUser_dl) >= 0];
end
F=[F,real(trace(sum(Q_dl,3))) <= Pd_W];
F=[F,q_ul >= 0];
F=[F,q_ul <= Pu_W];
F=[F,t_dl >= 1, t_ul >= 1]; % (24i)

%% Generate initial points
%% Iterative process        
obj=geomean([t_dl;t_ul]); % the objective in (24a)
MAX_ITER = 20; % max. number of iterations
lowerbound = zeros(MAX_ITER,1);
sumrate_true = zeros(MAX_ITER,1);
tic
for iIter = 1:MAX_ITER

    F1 = [];
    for iUser_dl=1:DL.NumUsers
        %/Constraint:
        RHS24b = EffNoisePower_dl; % the right-hand-side of (24b)
        RHS24b = RHS24b+real(DL.channel(:,iUser_dl)'*sum(Q_dl,3)*DL.channel(:,iUser_dl));
        RHS24b = RHS24b + sum(q_ul.*abs(g_ul_dl(:,iUser_dl)).^2);
        
        % Constraint (24b)
        F1=[F1,rcone([sqrt(1/(2*mypsi(iUser_dl)))*t_dl(iUser_dl);sqrt(mypsi(iUser_dl)/2)* mybeta(iUser_dl)],RHS24b,1/2)]; 

        % Constraint (24d)
        LHS24d = RHS24b - real(DL.channel(:,iUser_dl)'*Q_dl(:,:,iUser_dl)*DL.channel(:,iUser_dl)); % the left-hand-side of (24d)
        F1=[F1,LHS24d <= mybeta(iUser_dl)];

    end
    self_interference = EffNoisePower_ul*eye(Nr,Nr)+ Hli*sum(Q_dl,3)*Hli';

    for iUser_ul=1:UL.NumUsers
        Interf_ul = self_interference +UL.channel(:,iUser_ul+1:UL.NumUsers)*diag(q_ul(iUser_ul+1:UL.NumUsers))*UL.channel(:,iUser_ul+1:UL.NumUsers)';
        % Constraint (24c)
        F1=[F1,t_ul(iUser_ul)<= 1+ real(x_ul_init(iUser_ul)^2*UL.channel(:,iUser_ul)'*X_init(:,:,iUser_ul)^(-1) *UL.channel(:,iUser_ul))...
            + real((2*x_ul_init(iUser_ul)*UL.channel(:,iUser_ul)'*X_init(:,:,iUser_ul)^(-1) *UL.channel(:,iUser_ul))'*(x_ul(iUser_ul)-x_ul_init(iUser_ul)))...
            - real(trace((x_ul_init(iUser_ul)^2*(X_init(:,:,iUser_ul)^(-1))'*UL.channel(:,iUser_ul)*UL.channel(:,iUser_ul)'*(X_init(:,:,iUser_ul)^(-1))')'*...
            (Interf_ul-X_init(:,:,iUser_ul))))];

        %/Constraint (24e)
        F1=[F1,cone([2*x_ul(iUser_ul);1-q_ul(iUser_ul)], 1+q_ul(iUser_ul))];
    end

    diagnotics=solvesdp([F,F1],-obj,ops);
    if(diagnotics.problem==0) % successfully solved, optimal solution found
        lowerbound(iIter) = real((UL.NumUsers+DL.NumUsers)*log2(double(obj)));
        q_ul_init = double(q_ul);
        Q_dl_init = double(Q_dl);
        x_ul_init = double(x_ul);
        mypsi = double(t_dl./mybeta);
        self_interference = EffNoisePower_ul*eye(Nr,Nr)+ Hli*sum(Q_dl_init,3)*Hli';

        for iUser_ul=1:UL.NumUsers
            X_init(:,:,iUser_ul) = self_interference+UL.channel(:,iUser_ul+1:UL.NumUsers)...
                *diag(q_ul_init(iUser_ul+1:UL.NumUsers))*UL.channel(:,iUser_ul+1:UL.NumUsers)';
        end
        sumrate_true(iIter) = ComputeRates(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul);
    elseif ((3 <= diagnotics.problem) && (diagnotics.problem<=5))% some numerical issues, but could continue
        disp(strcat(diagnotics.info," at iteration ",num2str(iIter)))
        disp('Try using the current solution to continue')
        lowerbound(iIter) = real((UL.NumUsers+DL.NumUsers)*log2(double(obj)));
        q_ul_init = double(q_ul);
        Q_dl_init = double(Q_dl);
        x_ul_init = double(x_ul);
        mypsi = double(t_dl./mybeta);
        self_interference = EffNoisePower_ul*eye(Nr,Nr)+ Hli*sum(Q_dl_init,3)*Hli';

        for iUser_ul=1:UL.NumUsers
            X_init(:,:,iUser_ul) = self_interference+UL.channel(:,iUser_ul+1:UL.NumUsers)...
                *diag(q_ul_init(iUser_ul+1:UL.NumUsers))*UL.channel(:,iUser_ul+1:UL.NumUsers)';
        end
        sumrate_true(iIter) = ComputeRates(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul);
    else
        disp(strcat(diagnotics.info,' at iteration ',num2str(iIter)))
        sumrate_true(iIter:end) = [];
        lowerbound(iIter:end) = [];
        break
    end

end
toc
figure
plot(lowerbound)
hold on
plot(10:20,ones(11,1)*sumrate_true(20))
xlabel('Sensing SINR');
ylabel('Objective value');
ylim([0,100]);
legend('Lower bound of SumRate','True Sum Rate','Location','southeast')
% saveas(gcf, '../../results/ConvergencePlot_Algorithm2.png')
%%
result = ones(9,1)*sumrate_true(20);
sinr_required = transpose(7:1:15);
m = [ sinr_required result]
%%
DLSINR(1) = pow2db(DL_SINR(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul,1));
DLSINR(2) = pow2db(DL_SINR(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul,2));
ULSINR(1) = pow2db(UL_SINR(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul,1));
ULSINR(2) = pow2db(UL_SINR(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul,2));
%%
real(trace(sum(Q_dl_init,3)))
%%
% writematrix(m, '20241002_ComOnly_sinr_vs_capacity.csv');