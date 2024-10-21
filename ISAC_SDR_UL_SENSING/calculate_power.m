function power = calculate_power(VC_AMP, VS_PHASE) 
    power = norm(VC_AMP,"fro")^2 + norm(VS_PHASE,"fro")^2 ;
end