function OFDM = OFDM_Init()
    rng('default');                     % Set random number generator for reproducibility
    fc = 39e9;                          % Carrier frequency (Hz)
    c = 3e8;                            % Speed of light
    OFDM.N = 100;              % Number of subcarriers
    OFDM.delta_f = 0.1e9;       % Carrier spacing
    OFDM.carrier_freq = fc:OFDM.delta_f:fc+OFDM.N*OFDM.delta_f;
    OFDM.lambda = c./OFDM.carrier_freq;
end