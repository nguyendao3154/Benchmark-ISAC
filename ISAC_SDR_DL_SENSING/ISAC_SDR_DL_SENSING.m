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
 Pt = 50;                           % Peak transmit power (W)
 Nt =  32;                           % Number of array elements

sUCA = antenna_geometry(Nt, 39e9);
steervector = phased.SteeringVector('IncludeElementResponse', true, "SensorArray",sUCA);
%% OFDM params
OFDM = OFDM_Init();
%% Downlink
DL = DL_Init(steervector, fc, OFDM);
%% Radar
Radar = Radar_Init(steervector, fc, OFDM);
%%
sensitivity = 10; 
for gamma = 10:2:20
previous_result = 0;
VC_previous = rand(Nt,Nt,DL.NumUsers);
VS_previous = rand(Nt,Nt,Radar.NumTgts);
% VC_previous = a;
% VS_previous = b;

first_flag = 1;
stop_flag = 0;
count = 0;
iteration = 1;
% result(1) = compute_capacity(DL.channel, abs(v_tx), angle(v_tx));
% while stop_flag ~= 1
tic
% gamma = 10;
while iteration < 6
    cvx_begin sdp 
        % cvx_solver_settings('verbose', false); % or use specific solver settings
        % cvx_solver_settings('printlevel', 1);  % depends on solver, some solvers use 'printlevel' 
        % [r, res] = mosekopt('minimize echo(0)', prob);
        
        variable VC(Nt, Nt, DL.NumUsers) complex semidefinite  
        variable VS(Nt, Nt, Radar.NumTgts) complex semidefinite 
        % variable phase(Nt, Nt, DL.NumUsers)
        obj = 0;

        AC_previous = sum(VC_previous,3);
        AS_previous = sum(VS_previous,3);
        
        AC = sum(VC,3);
        AS = sum(VS,3);

        q_previous = (zeros(1, DL.NumUsers)) ;
        B = cvx(zeros(Nt,Nt, DL.NumUsers)) ;
        
        for j = 1:DL.NumUsers
            q_previous(j) = log(real(DL.channel(:,j)'*(sum(VC_previous(:,:,1:DL.NumUsers~=j),3) + sum(VS_previous,3))*DL.channel(:,j))+1)/log(2);
            B(:,:,j) = sum(VC(:,:,1:DL.NumUsers~=j),3) - sum(VC_previous(:,:,1:DL.NumUsers~=j),3) ;
        end


        for j = 1:DL.NumUsers
                obj = obj + log(real(DL.channel(:,j)'*(AC+AS)*DL.channel(:,j))+1)/log(2) ...
                   - q_previous(j) - ...
                   real(DL.channel(:,j)'*(B(:,:,j)+(AS-AS_previous))*DL.channel(:,j))/(2^(q_previous(j))*log(2));      
        end
        
        power = trace(sum(real(VC),3)) + trace(sum(real(VS),3));

        D_previous = (zeros(Nt,Nt, Radar.NumTgts));
        D = cvx(zeros(Nt,Nt, Radar.NumTgts));

        for m = 1:Radar.NumTgts
            D_previous(:,:,m) = Interference_Plus_Noise_Radar(Radar.combined, Radar.sv, VC_previous, VS_previous, m);
        end
        for m = 1:Radar.NumTgts
            D(:,:,m) = Interference_Plus_Noise_Radar(Radar.combined, Radar.sv, VC, VS, m);
        end
        maximize (obj)
        ratio = db2pow(gamma);
        subject to
            power <= 64;
            % trace(sum(real(VS),3)) >= 50;
            for j = 1:DL.NumUsers
                % real(DL.channel(:,j)'*((AC+AS)-(1+1/10)*VC(:,:,j))*DL.channel(:,j)) <=  1;
            end
            for m = 1:Radar.NumTgts
                real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*Radar.sv(:,m)) - real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*...
                (D(:,:,m) - D_previous(:,:,m))*inv(D_previous(:,:,m))*Radar.sv(:,m)) >= ...
                ratio/abs(Radar.combined(m))^2*real(inv_pos(Radar.sv(:,m)'*VS(:,:,m)*Radar.sv(:,m)));
                % real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*Radar.sv(:,m)) - real(Radar.sv(:,m)'*inv(D_previous(:,:,m))*...
                % (D(:,:,m) - D_previous(:,:,m))*inv(D_previous(:,:,m))*Radar.sv(:,m)) >= ...
                % 10/real(Radar.combined(m)^2)*real(inv(Radar.sv(:,m)'*VS(:,:,m)*Radar.sv(:,m)));
            end
    cvx_end
    if (strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed'))  ~= 1 
        VC_previous = VC;
        VS_previous = VS;
    else 
            VS_previous = rand(Nt,Nt,Radar.NumTgts);
            VC_previous = rand(Nt,Nt,DL.NumUsers);
            iteration = 1;
    end
    

    [eigvec, eigval] = eig(VC_previous(:,:,1));
    vc_tx(:,1) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

    [eigvec, eigval] = eig(VC_previous(:,:,2));
    vc_tx(:,2) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

    [eigvec, eigval] = eig(VS_previous(:,:,1));
    vs_tx(:,1) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

    [eigvec, eigval] = eig(VS_previous(:,:,2));
    vs_tx(:,2) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

    % [eigvec, eigval] = eig(VS_previous(:,:,3));
    % vs_tx(:,3) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);
    % 
    % [eigvec, eigval] = eig(VS_previous(:,:,4));
    % vs_tx(:,4) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);
    % 
    % result_SINR(1) = pow2db(DL_SINR(DL.channel, vc_tx, vs_tx, 1));
    % result_SINR(2) = pow2db(DL_SINR(DL.channel, vc_tx, vs_tx, 2));
    % result_SINR(3) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 1));
    % result_SINR(4) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 2));
    % result_capacity = (log(1+db2pow(result_SINR(1))) + log(1+db2pow(result_SINR(2))))/log(2);
    % 


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
    power_com(gamma+1) = trace(sum(real(VC),3));
    power_sen(gamma+1) = trace(sum(real(VS),3));
    result_SINR(1,gamma+1) = pow2db(DL_SINR(DL.channel, vc_tx, vs_tx, 1));
    result_SINR(2,gamma+1) = pow2db(DL_SINR(DL.channel, vc_tx, vs_tx, 2));

    result_SINR(3,gamma+1) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 1));
    result_SINR(4,gamma+1) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 2));
    % result_SINR(7,gamma+1) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 3));
    % result_SINR(8,gamma+1) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 4));
    best_result(gamma+1) = log(1+db2pow(result_SINR(1,gamma+1)))/log(2) + log(1+db2pow(result_SINR(2,gamma+1)))/log(2) ;
    % + log(1+db2pow(result_SINR(3,gamma+1)))/log(2) + log(1+db2pow(result_SINR(4,gamma+1)))/log(2);
    result_capacity(gamma+1) = cvx_optval;
end
toc
%%
plot(11:2:21, best_result(11:2:21));
ylim([10 35]);
a = VC_previous  ;
b = VS_previous  ;
matrix = [transpose(11:2:21) transpose(best_result(11:2:21))];
% writematrix(matrix, '20241002_DLSensing_sinr_capacity.csv');
%%
plot(0:20, best_result, 'LineWidth',2, 'Marker','o');
% ylim([0,100]);
    xlabel('Sensing SINR (dB)')
    ylabel('Total capacity (bit/s/Hz)')
%%
plot(0:10,result, 'LineWidth',2)
    % grid on
    xlabel('Iteration')
    ylabel('Objective value')
    set(gca,'FontSize', 12)
    ylim([0,50]);
    xlim([0,10]);
%%
plot(0:10,SINR(1,:),  'LineStyle','--', 'LineWidth',2, 'Marker','o'); hold on;
plot(0:10,SINR(2,:),  'LineStyle',':', 'LineWidth',2); hold on;
plot(0:10,SINR(3,:), 'LineStyle','-.', 'LineWidth',2, 'Marker','o'); hold on;
plot(0:10,SINR(4,:), 'LineStyle','-', 'LineWidth',2);
    xlabel('Iteration')
    ylabel('SINR (dB)')
xlim([0,10]);
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
%%
% Use the eig function to find eigenvalues and eigenvectors
[eigvec, eigval] = eig(VC_previous(:,:,1));
vc_tx(:,1) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

[eigvec, eigval] = eig(VC_previous(:,:,2));
vc_tx(:,2) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

[eigvec, eigval] = eig(VS_previous(:,:,1));
vs_tx(:,1) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

[eigvec, eigval] = eig(VS_previous(:,:,2));
vs_tx(:,2) = sqrt(eigval(Nt,Nt))*eigvec(:,Nt);

% v_tx = v_tx*v_tx';
%%
for j = 1:DL.NumUsers
    % vc_tx(:,j) = (DL.channel(:,j)'*VC(:,:,j)*DL.channel(:,j))^(-1/2)*VC(:,:,j)*DL.channel(:,j);
end

% VC_TX = 
total_pwr = norm(vc_tx,"fro")^2 + norm(vs_tx,"fro")^2 ;
pwr(1) = norm(vc_tx,"fro")^2 ;
pwr(2) = norm(vs_tx,"fro")^2 ;
result_SINR(1) = pow2db(DL_SINR(DL.channel, vc_tx, vs_tx, 1));
result_SINR(2) = pow2db(DL_SINR(DL.channel, vc_tx, vs_tx, 2));
result_SINR(3) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 1));
result_SINR(4) = pow2db(Radar_SINR(Radar.combined, Radar.sv, vc_tx, vs_tx, 2));
result_capacity = compute_capacity(DL.channel, vc_tx, vs_tx);
% result_capacity = compute_capacity(DL.channel, V);
%%
DL.channel(:,1)'*vc_tx(:,1)*vc_tx(:,1)'*DL.channel(:,1)
DL.channel(:,1)'*VC_previous(:,:,1)*DL.channel(:,1)
%%
result_SINR(1) = pow2db(DL_SINR(DL.channel, VC, VS, 1));
result_SINR(2) = pow2db(DL_SINR(DL.channel, VC, VS, 2));
result_SINR(3) = pow2db(Radar_SINR(Radar.combined, Radar.sv, VC, VS, 1));
result_SINR(4) = pow2db(Radar_SINR(Radar.combined, Radar.sv, VC, VS, 2));
result_capacity = compute_capacity(DL.channel, VC, VS);
