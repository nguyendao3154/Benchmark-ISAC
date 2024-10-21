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
 load('array_initialization.mat');
 Pt = 100;                           % Peak transmit power (W)
 Nt =  32;                           % Number of array elements

sUCA = antenna_geometry(Nt, 39e9);
steervector = phased.SteeringVector('IncludeElementResponse', true, "SensorArray",sUCA);
%% OFDM params
OFDM = OFDM_Init();
%% Uplink
UL = UL_Init(steervector, fc, OFDM);
%% Radar
Radar = Radar_Init(steervector, fc, OFDM);

%%
sensitivity = 10;
tic
for gamma = 10:2:20
previous_result = 0;
VS_previous = rand(Nt,Nt,Radar.NumTgts);
xc_previous = rand(UL.NumUsers,1);
e_previous = rand(UL.NumUsers,1);

power = trace(sum(real(VS_previous(:,:,1:Radar.NumTgts)),3));

first_flag = 1;
stop_flag = 0;
count = 0;


% result(1) = compute_capacity(DL.channel, abs(v_tx), angle(v_tx));
% while stop_flag ~= 1
iteration = 1;
% gamma = 15;
while iteration < 6
    cvx_begin sdp 
        % cvx_solver_settings('verbose', false); % or use specific solver settings
        % cvx_solver_settings('printlevel', 1);  % depends on solver, some solvers use 'printlevel' 
        % [r, res] = mosekopt('minimize echo(0)', prob);
        
        variable VS(Nt, Nt, Radar.NumTgts) complex semidefinite
        variable e(UL.NumUsers) nonnegative
        variable xc(UL.NumUsers) nonnegative
        variable u(UL.NumUsers) nonnegative
        % variable phase(Nt, Nt, DL.NumUsers)

        AS_previous = sum(VS_previous,3);
        
        AS = sum(VS,3);

        E_previous = (zeros(Nt,Nt, UL.NumUsers));
        E = cvx(zeros(Nt,Nt, UL.NumUsers));
        for k = 1:UL.NumUsers
            E_previous(:,:,k) = Interference_Plus_Noise_UL(Radar.combined, Radar.sv, VS_previous, UL.channel, e_previous, k);
        end
        for k = 1:UL.NumUsers
            E(:,:,k) = Interference_Plus_Noise_UL(Radar.combined, Radar.sv, VS, UL.channel, e, k);
        end

        obj = sum(log(real(1+u)));
        
        power = trace(sum(real(VS),3));

        D_previous = (zeros(Nt,Nt, Radar.NumTgts));
        D = cvx(zeros(Nt,Nt, Radar.NumTgts));

        for m = 1:Radar.NumTgts
            D_previous(:,:,m) = Interference_Plus_Noise_Radar(Radar.combined, Radar.sv, VS_previous, UL.channel, e_previous, m); 
        end
        for m = 1:Radar.NumTgts
            D(:,:,m) = Interference_Plus_Noise_Radar(Radar.combined, Radar.sv, VS, UL.channel, e, m);
        end

        maximize (obj)
        ratio = db2pow(gamma);
        subject to
            power <= 64;
            for k = 1:UL.NumUsers
                e(k) <= 3;
            end
            for m = 1:Radar.NumTgts
                real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*Radar.sv(:,m)) - real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*...
                (D(:,:,m) - D_previous(:,:,m))*inv(D_previous(:,:,m))*Radar.sv(:,m)) >= ...
                ratio/abs(Radar.combined(m))^2*real(inv_pos(Radar.sv(:,m)'*VS(:,:,m)*Radar.sv(:,m)));

                % real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*Radar.sv(:,m)) - real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*...
                % (D(:,:,m) - D_previous(:,:,m))*inv(D_previous(:,:,m))*Radar.sv(:,m)) >= ...
                % ratio/abs(Radar.combined(m))^2*real(inv(Radar.sv(:,m)'*VS(:,:,m)*Radar.sv(:,m)));
            end

            for k=1:UL.NumUsers
                % xc_previous(k)^2 + 2*xc_previous(k)*(xc(k) - xc_previous(k)) >= 10;
            end
            for k=1:UL.NumUsers
                u(k) <= xc_previous(k)^2 + 2*xc_previous(k)*(xc(k) - xc_previous(k));
            end

            for k =1:UL.NumUsers
                quad_over_lin(xc(k),e(k)) <= (real(UL.channel(:,k)'*inv(E_previous(:,:,k))*UL.channel(:,k)) - real(UL.channel(:,k)'*inv(E_previous(:,:,k))*...
                (E(:,:,k) - E_previous(:,:,k))*inv(E_previous(:,:,k))*UL.channel(:,k)));
            end
            
    cvx_end
    

    % 
    % result_SINR(1) = pow2db(UL_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 1));
    % result_SINR(2) = pow2db(UL_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 2));
    % result_SINR(3) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 1));
    % result_SINR(4) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 2));
    % result_capacity = cvx_optval;
    if strcmp(cvx_status,'Infeasible') ~= 1 
        VS_previous = VS;
        e_previous = e;
        xc_previous = xc;
    else 
            VS_previous = rand(Nt,Nt,Radar.NumTgts);
            xc_previous = rand(UL.NumUsers,1);
            e_previous = rand(UL.NumUsers,1);
            iteration = 1;
    end

    [eigvec, eigval] = eig(VS_previous(:,:,1));
    vs_tx(:,1) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);
    
    [eigvec, eigval] = eig(VS_previous(:,:,2));
    vs_tx(:,2) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);
    % sensitivity = abs(result_capacity - previous_result);
    % 
    % if sensitivity < 1e-2
    %     count = count + 1;
    %     if count == 5
    %         stop_flag = 1;
    %     end
    % else  
    %     count = 0;            
    % end
  
    % result(iteration) = result_capacity;
    % SINR(:,iteration) = result_SINR;
    % previous_result = result_capacity;
    iteration = iteration + 1;
end
    power_sen(gamma+1) = trace(sum(real(VS),3));
    result_SINR(1,gamma+1) = pow2db(UL_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 1));
    result_SINR(2,gamma+1) = pow2db(UL_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 2));
    result_SINR(3,gamma+1) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 1));
    result_SINR(4,gamma+1) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vs_tx, UL.channel, e_previous, 2));
    best_result(gamma+1) = log(1+db2pow(result_SINR(1,gamma+1)))/log(2) + log(1+db2pow(result_SINR(2,gamma+1)))/log(2);
end
toc
%%
plot(11:2:21, best_result(11:2:21));
% ylim([32,33]);
    xlabel('Sensing SINR')
    ylabel('Objective value')
matrix = [transpose(11:2:21) transpose(best_result(11:2:21))];
% writematrix(matrix, '20241002_ULSensing_sinr_capacity.csv');
%%
plot(0:5,result, 'LineWidth',2)
    % grid on
    xlabel('Iteration')
    ylabel('Objective value')
    set(gca,'FontSize', 12)
    ylim([0,50]);
    xlim([0,5]);
%%
plot(0:9,SINR(1,:),  'LineStyle','--', 'LineWidth',2, 'Marker','o'); hold on;
plot(0:9,SINR(2,:),  'LineStyle',':', 'LineWidth',2); hold on;
plot(0:9,SINR(3,:), 'LineStyle','-.', 'LineWidth',2, 'Marker','o'); hold on;
plot(0:9,SINR(4,:), 'LineStyle','-', 'LineWidth',2);
    xlabel('Iteration')
    ylabel('SINR (dB)')
xlim([0,9]);
    legend('User 1', 'User 2','Target 1', 'Target 2');
%%
Result_after_decompose = log(1+db2pow(SINR(1,:)))/log(2) + log(1+db2pow(SINR(2,:)))/log(2)
plot(0:9,Result_after_decompose, 'LineStyle','--', 'LineWidth',2); hold on;
% plot(0:6,result, 'LineWidth',2);
    % grid on
    xlabel('Iteration')
    ylabel('Objective value')
    set(gca,'FontSize', 12)
    ylim([0,50]);
    % xlim([0,5]);
    legend('After decompositon', 'Before decompositon');