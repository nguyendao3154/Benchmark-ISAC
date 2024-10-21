% function capacity = compute_capacity(H, VC, VS)
%     NumDL = size(VC, 3);
% 
%     dl_capacity = 0;
%     for j = 1:NumDL
%         dl_capacity = dl_capacity + log(1 + DL_SINR(H, VC, VS, j))/log(2);
%     end
% 
%     capacity = dl_capacity;
% end
%% Small V
function capacity = compute_capacity(H, VC, VS)
    NumDL = size(VC, 2);

    dl_capacity = 0;
    for j = 1:NumDL
        dl_capacity = dl_capacity + log(1 + DL_SINR(H, VC, VS, j))/log(2);
    end

    capacity = dl_capacity;
end
%% GA version
% function capacity = compute_capacity(H, VC_AMP, VC_PHASE, VS_AMP, VS_PHASE)
%     VC = VC_AMP.*(1j*VC_PHASE);
%     VS = VS_AMP.*(1j*VS_PHASE);
%     NumDL = size(VC, 2);
% 
%     dl_capacity = 0;
%     for j = 1:NumDL
%         dl_capacity = dl_capacity + log(1 + DL_SINR(H, VC_AMP, VC_PHASE, VS_AMP, VS_PHASE, j))/log(2);
%     end
% 
%     capacity = -dl_capacity;
% end

