%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file    ISAC_GreedySearch.m
% @author  Nguyen Dao - RS Group - UTwente.
% @version 1.0
% @date    Jan 09, 2024
% @brief   
% @history
% 
%                      Revision History                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revision     Date            By              Description                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.0.0        09-Jan-2024     Nguyen Dao          create                                %

%% Some constants
 Pt = 100;                           % Peak transmit power (W)
 Nt =  64;                           % Number of array elements
 load('array_initialization.mat');
%% Some more fixed parameters
%% OFDM params
OFDM = OFDM_Init();
%% Downlink
DL = DL_Init(steervector, fc, OFDM);
%% Radar
Radar = Radar_Init(steervector, fc, OFDM);

%% Calculate steering vector
ang_elv = linspace(0, 90, 181);       % Grid of Elevation angles
ang_az = linspace(-180, 180, 181);       % Grid of Elevation angles

sv_azimuth.DL1 = steervector(fc, [ang_az; DL.elevation(1).*ones(size(ang_az))]);
sv_elevation.DL1 = steervector(fc, [DL.azimuth(1).*ones(size(ang_elv)); ang_elv]);

sv_azimuth.DL2 = steervector(fc, [ang_az; DL.elevation(2).*ones(size(ang_az))]);
sv_elevation.DL2 = steervector(fc, [DL.azimuth(2).*ones(size(ang_elv)); ang_elv]);

%%
for i = 1:181
    sv_3d(:,i,:) = steervector(fc, [ang_az(i)*ones(1,181); ang_elv]);
end
%% Spectral Efficiency Optimization
H_combined_all = [DL.channel];
precoder = H_combined_all*inv(H_combined_all'*H_combined_all);

power = sqrt(norm(precoder,"fro"))^2;
v_tx(:,1) = precoder(:,1)/power*sqrt(100);
v_tx(:,2) = precoder(:,2)/power*sqrt(100);

new_pwr = norm(v_tx,"fro")^2 ;

for j = 1:DL.NumUsers
    V_previous(:,:,j) = (v_tx(:,j)*v_tx(:,j)');
end

power_0 = 0;
for j = 1:DL.NumUsers
            power_0 = power_0 + trace(V_previous(:,:,j));
end
%%
tic
vc_amp = optimvar('vc_amp', Nt, DL.NumUsers, 'LowerBound',0, 'UpperBound', 5);
vc_phase = optimvar('vc_phase', Nt, DL.NumUsers, 'LowerBound',0, 'UpperBound', 2*pi);

vs_amp = optimvar('vs_amp', Nt, Radar.NumTgts, 'LowerBound',0, 'UpperBound', 5);
vs_phase = optimvar('vs_phase', Nt, Radar.NumTgts, 'LowerBound',0, 'UpperBound', 2*pi);

obj = fcn2optimexpr( @compute_capacity, DL.channel, vc_amp, vc_phase, vs_amp, vs_phase);

% parfor gamma = 11:12

prob = optimproblem('Objective',obj);

prob.Constraints.confn2_1 = ...
    fcn2optimexpr( @DL_SINR, DL.channel, vc_amp, vc_phase, vs_amp, vs_phase, 1) >= db2pow(10);

prob.Constraints.confn2_2 = ...
    fcn2optimexpr( @DL_SINR, DL.channel, vc_amp, vc_phase, vs_amp, vs_phase, 2) >= db2pow(10);

prob.Constraints.confn3_1 = ...
    fcn2optimexpr( @Radar_SINR, Radar.combined, Radar.sv, vc_amp, vc_phase, vs_amp, vs_phase, 1) >= db2pow(10);

prob.Constraints.confn3_2 = ...
    fcn2optimexpr( @Radar_SINR, Radar.combined, Radar.sv, vc_amp, vc_phase, vs_amp, vs_phase, 2) >= db2pow(10);


prob.Constraints.confn4 = ...
    fcn2optimexpr(@calculate_power, vc_amp, vs_amp) <= 100;

options = optimoptions('ga','PlotFcn', {@gaplotbestf, @gaplotscores});
options.PopulationSize = 200;
options.EliteCount = 10;
options.ConstraintTolerance = 0;
options.FunctionTolerance = 0;
options.MaxGenerations = 20;
options.MaxStallGenerations = 50;
options.SelectionFcn = {@selectiontournament,2};
options.UseParallel =  true;
% options.CrossoverFcn = {@crossoverscattered};
% options.MutationFcn = {@mutationadaptfeasible};
[fit2,fval] = solve(prob,'Solver','ga','Options',options);
toc
% end
%%
plot(-brushedData(:,2),'LineWidth',2);
    % grid on
    xlabel('Iteration')
    ylabel('Objective value')
    set(gca,'FontSize', 12)
%%
result_SINR(1) = pow2db(DL_SINR(DL.channel, fit2.vc_amp, fit2.vc_phase, fit2.vs_amp, fit2.vs_phase, 1));
result_SINR(2) = pow2db(DL_SINR(DL.channel, fit2.vc_amp, fit2.vc_phase, fit2.vs_amp, fit2.vs_phase, 2));
result_SINR(3) = pow2db(Radar_SINR(Radar.combined, Radar.sv, fit2.vc_amp, fit2.vc_phase, fit2.vs_amp, fit2.vs_phase, 1));
result_SINR(4) = pow2db(Radar_SINR(Radar.combined, Radar.sv, fit2.vc_amp, fit2.vc_phase, fit2.vs_amp, fit2.vs_phase, 2));
result_capacity = compute_capacity(DL.channel, fit2.vc_amp, fit2.vc_phase, fit2.vs_amp, fit2.vs_phase);
% new_pwr = calculate_power(v_tx) ;
   %%
   figure;
   % v_tx(:,1) = v_tx(:,1)/norm(v_tx(:,1));
    for i = 1:181
        Pattern_elv1(i) = abs(sv_elevation.DL2(:,i)'*(v_tx(:,2)/8))^2;
    end
    hold on
    grid on
    xlabel('Elevation (deg)')
    ylabel('(dB)')
    xlim([0 90])
    ylim([-20 20])
    plot(ang_elv, pow2db(Pattern_elv1), 'LineWidth', 2, 'LineStyle','-.', 'Color','black')
    set(gca,'FontSize', 13)
    % legend('DL 1', 'DL 2', 'Sensing 1', 'Sensing 2', 'Sum TX', 'Location', 'southoutside', 'Orientation', 'horizontal')
    % title('Beam Pattern')
    %%
    % csvwrite('elv_patt.csv',transpose(pow2db(Pattern_elv_ul1)));
%% Desired DL pattern Elevation
    for i = 1:181
        Pattern_az1(i) = abs(sv_azimuth.DL1(:,i)'*v_tx(:,1))^2;
    end

    figure
    hold on
    grid on
    xlabel('Azimuth (deg)')
    ylabel('(dB)')
    xlim([-180 180])

    plot(ang_az, pow2db(Pattern_az1), 'LineWidth', 2, 'LineStyle','-.','Color','black')

    set(gca,'FontSize', 13)
    % legend('DL 1', 'DL 2', 'Sensing 1', 'Sensing 2', 'Sum TX', 'Location', 'southoutside', 'Orientation', 'horizontal')
    % title('Beam Pattern')
    %%
    % csvwrite('az_patt.csv',transpose(pow2db(Pattern_az_sum_rx)));
