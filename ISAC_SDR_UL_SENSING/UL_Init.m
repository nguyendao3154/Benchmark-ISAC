function UL = UL_Init(steervector, fc, OFDM)


% Uplink
UL.NumUsers = 2;                              % Number of communication users
UL.range = 200*ones(1, UL.NumUsers);
UL.azimuth = [0, 120];
UL.elevation = [60, 60]; 
UL.UsersPos = UL.range.*[cosd(UL.elevation).*sind(UL.azimuth); cosd(UL.elevation).*cosd(UL.azimuth); sind(UL.elevation)];
UL.sv = steervector(fc, [UL.azimuth; UL.elevation]);
NumElements = size(UL.sv,1);
UL.sv = UL.sv/sqrt(NumElements);
% Channel

    UL.Pathloss = (OFDM.lambda(1)/(4*pi)./UL.range(1))^2;

    % UL.channel = sqrt(NumElements)*OFDM.lambda(1)/(4*pi)./UL.range.*exp(-1j*2*pi.*UL.range./OFDM.lambda(1)).*UL.sv;
    % UL.channel = OFDM.lambda(1)/(4*pi)./UL.range.*exp(-1j*2*pi.*UL.range./OFDM.lambda(1)).*UL.sv;
    UL.channel = ones(1,UL.NumUsers)*sqrt(1e-12).*UL.sv/sqrt(1e-13);
end