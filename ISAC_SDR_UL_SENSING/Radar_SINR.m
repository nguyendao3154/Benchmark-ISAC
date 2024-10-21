%% SMALL V version
function radar_SINR = Radar_SINR(BETA, STEERING, vs, G, p, M) 
    NumTgts = size(vs, 2);

    for m = 1:NumTgts   
        VS(:,:,m) = vs(:,m)*vs(:,m)'; 
    end


    D = Interference_Plus_Noise_Radar(BETA, STEERING, VS, G, p, M);

    radar_SINR =  abs((BETA(M)*STEERING(:,M)'*vs(:,M))^2*(STEERING(:,M)') ...
        *inv(D)*(STEERING(:,M)));
end

%% BIG V version
% function radar_SINR = Radar_SINR(BETA, STEERING, VS, G, p, M) 
%     NumElemt = size(VS,1);
%     NumTgts = size(VS, 3);
% 
%     D = Interference_Plus_Noise_Radar(BETA, STEERING, VS, G, p, M);
% 
%     radar_SINR =  abs((BETA(M)^2*STEERING(:,M)'*VS(:,:,M))*STEERING(:,M)*(STEERING(:,M)') ...
%         *inv(D)*(STEERING(:,M)));
% end
%% GA version
% function radar_SINR = Radar_SINR(BETA, STEERING, VC_AMP, VC_PHASE, VS_AMP, VS_PHASE, M) 
%     VC = VC_AMP.*(1j*VC_PHASE);
%     VS = VS_AMP.*(1j*VS_PHASE);
%     NumElemt = size(VC,1);
%     NumTgts = size(VS, 2);
%     NumDL = size(VC, 2);
% 
%     D = Interference_Plus_Noise_Radar(BETA, STEERING, VC*VC', VS*VS', M);
% 
%     radar_SINR =  abs((BETA(M)*STEERING(:,M)'*VS(:,M))^2*(STEERING(:,M)') ...
%         *inv(D)*(STEERING(:,M)));
% end