function D = Interference_Plus_Noise_Radar(BETA, STEERING, VC, VS, M)
    NumElemt = size(VC,1);
    NumTgts = size(VS, 3);
    NumDL = size(VC, 3);

    DL_Interference = zeros(NumElemt,NumElemt);
    for m = 1:NumTgts
         DL_Interference = DL_Interference + (BETA(m)^2*STEERING(:,m)*STEERING(:,m)'*sum(VC(:,:,1:NumDL),3) ...
             *STEERING(:,m)*STEERING(:,m)');
    end

    Radar_Interference  = (BETA(M)^2*STEERING(:,M)*STEERING(:,M)'*sum(VS(:,:,1:NumTgts~=M),3)*STEERING(:,M)*STEERING(:,M)');
    for m_prime = 1:NumTgts
        if m_prime ~= M
            Radar_Interference = Radar_Interference + BETA(m_prime)^2*STEERING(:,m_prime)*STEERING(:,m_prime)'*...
                sum(VS(:,:,1:NumTgts),3)*STEERING(:,m_prime)*STEERING(:,m_prime)';
        end
    end

    Noise_power = eye(NumElemt)*1;

    D = Radar_Interference + DL_Interference + Noise_power;
end