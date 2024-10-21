function E = Interference_Plus_Noise_UL(BETA, STEERING, VS, G, p, K)
    NumElemt = size(VS,1);
    NumTgts = size(VS, 3);
    NumUL = size(p,1);
    
    UL_Interference = zeros(NumElemt,NumElemt);
    
    for k_prime = 1:NumUL
        if k_prime ~= K
            UL_Interference = UL_Interference + G(:,k_prime)*G(:,k_prime)'*p(k_prime);
        end
    end

    Radar_Interference = zeros(NumElemt,NumElemt);

    for m = 1:NumTgts
        Radar_Interference = Radar_Interference + BETA(m)^2*STEERING(:,m)*STEERING(:,m)'*sum(VS,3)*STEERING(:,m)*STEERING(:,m)';
    end

    Noise_power = eye(NumElemt)*1;

    E = Radar_Interference + UL_Interference + Noise_power;

end