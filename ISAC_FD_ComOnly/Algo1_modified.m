%*********************%
%  SYSTEM PARAMETERS  %
%*********************%
%% Some constants
rng('default');                     % Set random number generator for reproducibility
fc = 39e9;                          % Carrier frequency (Hz)
c = 3e8;                            % Speed of light
Pt = 200;                         % Peak transmit power (W)
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
DL.range = 100*ones(1, DL.NumUsers);
% DL.azimuth = [0];
% DL.elevation = [70];
DL.azimuth = [0, 120];
DL.elevation = [80, 80];
DL.UsersPos = DL.range.*[cosd(DL.elevation).*sind(DL.azimuth); cosd(DL.elevation).*cosd(DL.azimuth); sind(DL.elevation)];
DL.sv = steervector(fc, [DL.azimuth; DL.elevation]);

% Channel
AntennaGain_default = 1;
    DL.Pathloss = (OFDM.lambda(1)/(4*pi)*...
        AntennaGain_default/DL.range(1))^2;
    DL.channel = OFDM.lambda(1)/(4*pi)*...
        AntennaGain_default./DL.range.*exp(-1j*2*pi.*DL.range./OFDM.lambda(1)).*DL.sv;
    %% Uplink
UL.NumUsers = 2;                              % Number of communication users
UL.range = 100*ones(1, DL.NumUsers);
% UL.azimuth = [0];
% UL.elevation = [50];
UL.azimuth = [0, 120];
UL.elevation = [60, 60];
UL.UsersPos = DL.range.*[cosd(UL.elevation).*sind(UL.azimuth); cosd(UL.elevation).*cosd(UL.azimuth); sind(UL.elevation)];
UL.sv = steervector(fc, [UL.azimuth; UL.elevation]);

% Channel
    UL.Pathloss = OFDM.lambda(1)/(4*pi)*...
        AntennaGain_default/UL.range(1);
    UL.channel = OFDM.lambda(1)/(4*pi)*...
        AntennaGain_default./UL.range.*exp(-1j*2*pi.*UL.range./OFDM.lambda(1)).*UL.sv;
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
NoisePower_dl = 1e-13;
NoisePower_ul = 1e-13;
scale_dl = NoisePower_dl;
DL.channel = DL.channel/sqrt(scale_dl); % scale DL channels with scaling factor of scale_dl
EffNoisePower_dl = NoisePower_dl/scale_dl; % the effective noise power after DL channels scaling

scale_ul = NoisePower_ul; 
UL.channel = UL.channel/sqrt(scale_ul); % scale UL channels with scaling factor of scale_ul 
Hli = Hli/sqrt(scale_ul); % scale SI channel with scaling factor of scale_ul
EffNoisePower_ul = NoisePower_ul/scale_ul; % effective noise power after UL channels scaling

%% Iterative MAXDET algorithm
% generate initial points
Pu_W = 50;
Pd_W = 200;
q_ul_init = Pu_W;

g_ul_dl = zeros(UL.NumUsers, DL.NumUsers);

Q_dl_init_tmp = zeros(Nt,Nt,DL.NumUsers);     % beamformer
for iUser_dl=1:DL.NumUsers
    a = rand(Nt,Nt)+1i*rand(Nt,Nt);
    Q_dl_init_tmp(:,:,iUser_dl) = a*a';
end
P_tmp = real(trace(sum(Q_dl_init_tmp,3)));
Q_dl_init = Q_dl_init_tmp/(P_tmp/(Pd_W));

Q_dl = sdpvar(Nt,Nt,DL.NumUsers,'hermitian','complex'); % DL covariance matrices variables
q_ul = sdpvar(UL.NumUsers,1); % UL power variables
ops=sdpsettings('solver','sdpt3','verbose',0); % solver options
% define constraints:
F=[];
for iUser_dl=1:DL.NumUsers
    F=[F,Q_dl(:,:,iUser_dl) >= 0];  % positive semidefinite matrix constraint
end
F=[F,real(trace(sum(Q_dl,3))) <= Pd_W]; % maximum power total constraint
F=[F,q_ul >= 0];        
F=[F,q_ul <= Pu_W];         % power constraint each user
MAX_ITER = 30; % max. number of iterations
lowerbound = zeros(MAX_ITER,1);
sumrate_true = zeros(MAX_ITER,1);
for iIter = 1:MAX_ITER
    sumrate_approximate = 0;
    % DL rates approximation
    for iUser_dl=1:DL.NumUsers
        Num_dl = EffNoisePower_dl+(DL.channel(:,iUser_dl)'*sum(Q_dl,3)*DL.channel(:,iUser_dl))...
                        +sum(q_ul.*(abs(g_ul_dl(:,iUser_dl))).^2);

        Den_dl_const = EffNoisePower_dl+real(DL.channel(:,iUser_dl)'*sum(Q_dl_init(:,:,1:DL.NumUsers~=iUser_dl),3)*DL.channel(:,iUser_dl))...
                        +sum(q_ul_init.*(abs(g_ul_dl(:,iUser_dl))).^2);
        
        Den_dl = real(((1/Den_dl_const)*DL.channel(:,iUser_dl)'...
            *sum(Q_dl(:,:,1:DL.NumUsers~=iUser_dl)-Q_dl_init(:,:,1:DL.NumUsers~=iUser_dl),3))...
            *DL.channel(:,iUser_dl))+(1/Den_dl_const)*sum((q_ul-q_ul_init).*(abs(g_ul_dl(:,iUser_dl))).^2);

        sumrate_approximate = sumrate_approximate + (logdet(Num_dl)) - real(log(Den_dl_const)) - Den_dl;
    end

    % UL rates approximation
    Num_ul = EffNoisePower_ul*eye(Nr)+Hli*sum(Q_dl,3)*Hli';
    Num_ul = Num_ul + UL.channel*diag(q_ul)*UL.channel';
    Den_ul_const = EffNoisePower_ul*eye(Nr) + Hli*sum(Q_dl_init,3)*Hli';
    
    Den_ul = real(trace((Hli'*Den_ul_const^(-1)*Hli)*sum(Q_dl-Q_dl_init,3)));
    
    sumrate_approximate = sumrate_approximate + logdet(Num_ul)- real(log(det(Den_ul_const))) - Den_ul;


    diagnotics=solvesdp(F,-sumrate_approximate,ops);
    if(diagnotics.problem==0) % successfully solved, optimal solution found
        lowerbound(iIter) = real(double(sumrate_approximate*log2(exp(1))));
        q_ul_init = double(q_ul);
        Q_dl_init = double(Q_dl);
        sumrate_true(iIter) = ComputeRates(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul);
    elseif ((3 <= diagnotics.problem) && (diagnotics.problem<=5)) % some numerical issues, but could continue
        disp(strcat(diagnotics.info," at iteration ",num2str(iIter)))
        disp('Try using the current solution to continue')
        lowerbound(iIter) = real(double(sumrate_approximate*log2(exp(1))));
        q_ul_init = double(q_ul);
        Q_dl_init = double(Q_dl);
        sumrate_true(iIter) = ComputeRates(DL.channel,UL.channel,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul);
    else
        sumrate_true(iIter:end) = [];
        lowerbound(iIter:end) = [];
        break
    end
            
end
figure
plot(lowerbound)
hold on
plot(sumrate_true)
legend('Lower bound of SumRate','True Sum Rate','Location','southeast')
% saveas(gcf, '../../results/ConvergencePlot_Algorithm1.png')