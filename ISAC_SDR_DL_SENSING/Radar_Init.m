function Radar = Radar_Init(steervector, fc, OFDM)   
    c = 3e8;
    Radar.NumTgts = 2;                      % Number of targets
    % Radar.azimuth = [0];
    % Radar.elevation = [30];
    Radar.azimuth = [0, 120];
    Radar.elevation = [70, 70];
    % Radar.azimuth = [0, 60, 120, 180];
    % Radar.elevation = [60, 60, 60, 60];
    Radar.sv = steervector(fc, [Radar.azimuth; Radar.elevation]);
    NumElements = size(Radar.sv,1);
    Radar.sv = Radar.sv/sqrt(NumElements);

    % Radar.pathloss = (4*pi)^(3/2)*(Radar.range).^2*fc/c.*exp(-1j*2*2*pi.*Radar.range./OFDM.lambda(1));
    Radar.cross_section = 1*ones(1,Radar.NumTgts);
    Radar.combined = ones(1,Radar.NumTgts)*sqrt(5);

    % Radar.combined = Radar.cross_section./Radar.pathloss;
    % Radar.combined = ones(1, Radar.NumTgts)*1e-8;
end