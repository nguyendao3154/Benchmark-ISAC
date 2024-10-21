%% Small V version
function DL_SINR = DL_SINR(H, vc, vs, J)
    NumDL = size(vc,2);
    NumTgt = size(vs, 2);
    DL_MUI = 0;

    for j = 1:NumDL
        if j ~= J
            DL_MUI = DL_MUI + real(H(:,J)'*vc(:,j)*vc(:,j)'*H(:,J));
        end
    end

    interference = 0;

    for m = 1:NumTgt
        interference = interference + real(H(:,J)'*vs(:,m)*vs(:,m)'*H(:,J));
    end
    DL_SINR = real(H(:,J)'*vc(:,J)*vc(:,J)'*H(:,J))/(DL_MUI + interference +1);
end

%% BIG V version
% function DL_SINR = DL_SINR(H, VC, VS, J)
%     NumDL = size(VC,3);
%     NumTgt = size(VS, 3);
%     DL_MUI = 0;
% 
%     for j = 1:NumDL
%         if j ~= J
%             DL_MUI = DL_MUI + real(H(:,J)'*VC(:,:,j)*H(:,J));
%         end
%     end
%     interference = 0;
%     for m = 1:NumTgt
%         interference = interference + real(H(:,J)'*VS(:,:,m)*H(:,J));
%     end
%     DL_SINR = real(H(:,J)'*VC(:,:,J)*H(:,J))/(DL_MUI + interference +1e-13);
% end

%% GA version
% function DL_SINR = DL_SINR(H, VC_AMP, VC_PHASE, VS_AMP, VS_PHASE, J)
%     VC = VC_AMP.*(1j*VC_PHASE);
%     VS = VS_AMP.*(1j*VS_PHASE);
%     NumDL = size(VC,2);
%     NumTgt = size(VS, 2);
%     DL_MUI = 0;
% 
%     for j = 1:NumDL
%         if j ~= J
%             DL_MUI = DL_MUI + abs(H(:,J)'*VC(:,j)*VC(:,j)'*H(:,J));
%         end
%     end
%     interference = 0;
%     for m = 1:NumTgt
%         interference = interference + abs(H(:,J)'*VS(:,m)*VS(:,j)'*H(:,J));
%     end
%     DL_SINR = abs(H(:,J)'*VC(:,J)*VC(:,j)'*H(:,J))/(DL_MUI + interference +1e-13);
% end