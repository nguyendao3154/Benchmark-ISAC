function DL = DL_Init(steervector, fc, OFDM)

% Downlink
DL.NumUsers = 2;                              % Number of communication users
DL.azimuth = [0, 120];
DL.elevation = [80, 80]; 
DL.range = 500;
% DL.azimuth = [0, 60, 0,  90];
% DL.elevation = [80, 80, 60, 60]; 
% DL.UsersPos = DL.range.*[cosd(DL.elevation).*sind(DL.azimuth); cosd(DL.elevation).*cosd(DL.azimuth); sind(DL.elevation)];
DL.sv = steervector(fc, [DL.azimuth; DL.elevation]);
NumElements = size(DL.sv,1);
DL.sv = DL.sv/sqrt(NumElements);
% Channel

    % DL.Pathloss = ones(1,DL.NumUsers)*1e-13;
    DL.channel = ones(1,DL.NumUsers)*sqrt(1e-12).*DL.sv/sqrt(1e-13);
    % DL.channel = OFDM.lambda(1)/(4*pi)./DL.range.*exp(-1j*2*pi.*DL.range./OFDM.lambda(1)).*DL.sv./sqrt(5e-13);
end