
function ULSINR = UL_SINR(H_dl,H_ul,Hli,g_ul_dl,Q_dl,q_ul,NoisePower_dl,NoisePower_ul,K)
                        


nUser_ul = size(H_ul,2);      % The number of UL users

nRx = size(Hli,1);           % number of receive antennas at BS on UL

self_interference = Hli*sum(Q_dl,3)*Hli';

    Interf_ul = self_interference+NoisePower_ul*eye(nRx,nRx);
    Interf_ul = Interf_ul + H_ul(:,K+1:nUser_ul)*diag(q_ul(K+1:nUser_ul))*H_ul(:,K+1:nUser_ul)';

    % compute SINR for each user
    ULSINR=q_ul(K)*real(H_ul(:,K)'*Interf_ul^(-1)*H_ul(:,K));

            

