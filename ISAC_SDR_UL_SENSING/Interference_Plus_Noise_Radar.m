function D = Interference_Plus_Noise_Radar(BETA, STEERING, VS, G, p, M)
    NumElemt = size(VS,1);
    NumTgts = size(VS, 3);
    NumUL = size(p,1);
    
    UL_Interference = zeros(NumElemt,NumElemt);
    for k = 1:NumUL
        UL_Interference = UL_Interference + (G(:,k)*G(:,k)'*p(k));
    end

    Radar_Interference  = (BETA(M)^2*STEERING(:,M)*STEERING(:,M)'*sum(VS(:,:,1:NumTgts~=M),3)*STEERING(:,M)*STEERING(:,M)');

    for m_prime = 1:NumTgts
        if m_prime ~= M
            Radar_Interference = Radar_Interference + BETA(m_prime)^2*STEERING(:,m_prime)*STEERING(:,m_prime)'*...
                sum(VS(:,:,1:NumTgts),3)*STEERING(:,m_prime)*STEERING(:,m_prime)';
        end
    end

    Noise_power = eye(NumElemt)*1;

    D = Radar_Interference + UL_Interference + Noise_power;

end