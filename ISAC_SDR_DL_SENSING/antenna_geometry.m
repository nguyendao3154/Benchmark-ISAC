%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file    ISAC_GreedySearch.m
% @author  Nguyen Dao - RS Group - UTwente.
% @version 1.0
% @date    Jan 09, 2024
% @brief   
% @history
% 
%                      Revision History                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revision     Date            By              Description                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.0.0        09-Jan-2024     Nguyen Dao          create                                %

function ant_array = antenna_geometry(N, F)
    %% Some constants
    rng('default');                     % Set random number generator for reproducibility
    lambda = freq2wavelen(F);          % Wavelength (m)
    Pt = 0.5e6;                         % Peak transmit power (W)
    c = 3e8;                            % Speed of light
    %%
    azang = -180:180;
    elang = -90:90;
    Ele_pattern = [zeros(1,90) ones(1,91)];
    magpattern = mag2db(repmat(Ele_pattern',1,numel(azang)));
    phasepattern = zeros(size(magpattern));
    custom_element = phased.CustomAntennaElement('AzimuthAngles',azang, ...
        'ElevationAngles',elang,'MagnitudePattern',magpattern, ...
        'PhasePattern',phasepattern);
    %%
    azang = (0:N-1)*360/N-180;
    array_radius = lambda/(4*sind(180/N));
    azang = (0:N-1)*360/N-180;
    
    ant_array = phased.ConformalArray(...
        'ElementPosition',[array_radius*cosd(azang);array_radius*sind(azang);[0:(N-1)].*zeros(1,N)*lambda/8],...
        'ElementNormal',[azang;zeros(1,N)],'Element',custom_element);
    
        % ant_array = phased.URA('Size',[sqrt(N), sqrt(N)], 'ElementSpacing',[lambda/2, lambda/2], ...
        % 'ArrayNormal','z');

    viewArray(ant_array,'ShowIndex',[1 N],'ShowNormals',true, ...
        'ShowLocalCoordinates',true,'Orientation',[0;0;0], ...
        'ShowAnnotation',true)
    %%
    figure;
    pattern(ant_array,F,-180:180,0,'PropagationSpeed',c,...
        'CoordinateSystem','polar',...
        'Type','efield','Normalize',true)
    %%
    figure;
    pattern(ant_array,F,0,-90:90,'PropagationSpeed',c,...
        'CoordinateSystem','polar',...
        'Type','efield','Normalize',true)
    %%
    figure;
    pattern(ant_array,F,-180:180,-90:90,'PropagationSpeed',c,...
        'CoordinateSystem','polar',...
        'Type','efield','Normalize',true)
    
end