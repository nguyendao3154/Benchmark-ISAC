
function DLSINR = DL_SINR(H_dl,H_ul,Hli,g_ul_dl,Q_dl,q_ul,NoisePower_dl,NoisePower_ul,L)
                        

nUser_dl = size(H_dl,2);      % The number of DL users



    Interf_ul_dl = sum(q_ul.*abs(g_ul_dl(:,L)).^2);
    Interf_dl = NoisePower_dl + real(H_dl(:,L)'*Q_dl(:,:,1:nUser_dl~=L)*H_dl(:,L));
    DLSINR = real(H_dl(:,L)'*Q_dl(:,:,L)*H_dl(:,L))/(Interf_dl+Interf_ul_dl);


            

