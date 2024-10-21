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
    
    % %% Desired pattern Azimuth
    % tgtAz = [-45 150 60 -120];                % Azimuths of the targets of interest
    % tgtRng = [5.31e3, 6.23e3, 5.7e3];   % Ranges of the targets of interest
    % ang = linspace(-180, 180, 181);       % Grid of azimuth angles
    % beamwidth = 10;                     % Desired beamwidth
    % % Desired beam pattern
    % idx = false(size(ang));
    % for i = 1:numel(tgtAz)
    %     idx = idx | ang >= tgtAz(i)-beamwidth/2 & ang <= tgtAz(i)+beamwidth/2;
    % end
    % Bdes = zeros(size(ang));
    % Bdes(idx) = 1;
    % 
    % figure;
    % polarplot(ang/180*pi, Bdes, 'LineWidth', 2)
    % thetaticks([0:30:360])
    % thetaticklabels({'0','30','60','90','120','150','180','-150','-120','-90','-60','-30'})
    % % xlabel('Azimuth (deg)')
    % title('Desired Beam Pattern')
    % grid on
    % 
    % %%
    % % rxpos = sCA.getElementPosition();   
    % % normalizedPos = rxpos/lambda;
    % steervector = phased.SteeringVector("SensorArray",ant_array,'IncludeElementResponse', true)
    % sv = steervector(F, [ang; 10*ones(size(ang))]);
    % %%
    % cvx_begin
    %    variable R(N,N);
    %    obj = 0;
    %     for i = 1:181
    %         obj = obj + pow_abs(Bdes(i) - sv(:,i)'*R*sv(:,i),2);
    %     end
    %     % a = sum(obj);
    % 
    %     minimize(obj)
    %     subject to 
    % 
    %        R >= 0;
    % cvx_end
    % %%
    % obj_res = 0;
    % for i = 1:181
    %     Bmmse(i) = abs(sv(:,i)'*R*sv(:,i));
    %     obj_res = obj_res + pow_abs(Bdes(i) - sv(:,i)'*R*sv(:,i),2);
    % end
    % 
    % 
    % %%
    % figure
    % hold on
    % plot(ang, pow2db(Bdes + eps), 'LineWidth', 2)
    % plot(ang, pow2db(Bmmse/max(Bmmse)), 'LineWidth', 2)
    % 
    % grid on
    % xlabel('Azimuth (deg)')
    % ylabel('(dB)')
    % legend('Desired', 'CVX Covariance', 'Location', 'southoutside', 'Orientation', 'horizontal')
    % ylim([-20 1])
    % % xlim([-90 90])
    % title('Transmit Beam Pattern')
    % %% Desired pattern Elevation
    % tgtElv = [50 10];                % Azimuths of the targets of interest
    % tgtRng = [5.31e3, 6.23e3, 5.7e3];   % Ranges of the targets of interest
    % ang = linspace(0, 90, 181);       % Grid of Elevation angles
    % beamwidth = 10;                     % Desired beamwidth
    % % Desired beam pattern
    % idx = false(size(ang));
    % for i = 1:numel(tgtElv)
    %     idx = idx | ang >= tgtElv(i)-beamwidth/2 & ang <= tgtElv(i)+beamwidth/2;
    % end
    % Bdes = zeros(size(ang));
    % Bdes(idx) = 1;
    % 
    % figure;
    % polarplot(ang/180*pi, Bdes, 'LineWidth', 2)
    % thetaticks([0:30:360])
    % thetaticklabels({'0','30','60','90','120','150','180','-150','-120','-90','-60','-30'})
    % % xlabel('Azimuth (deg)')
    % title('Desired Beam Pattern')
    % grid on
    % 
    % %%  
    % steervector = phased.SteeringVector('IncludeElementResponse', true, "SensorArray",ant_array)
    % sv = steervector(F, [45*ones(size(ang)); ang]);
    % 
    % %%
    % 
    % cvx_begin
    %     variable R(N,N)
    %     obj = 0;
    %     for i = 1:181
    %         obj = obj + pow_abs(Bdes(i) - sv(:,i)'*R*sv(:,i),2);
    %     end
    % 
    %     minimize( obj )
    %     subject to 
    %        for k = 1:N-1
    %             R(k+1,k+1) == R(k,k) 
    %        end
    % 
    %        % for i = 1:91
    %        %  abs(sv(:,i)'*R*sv(:,i) - Bdes(i)) <= 0.5
    %        % end
    %        R >= 0;
    % cvx_end
    % 
    % %%
    % for i = 1:181
    % Bmmse(i) = abs(sv(:,i)'*R*sv(:,i));
    % end
    % 
    % %%
    % figure
    % hold on
    % plot(ang, pow2db(Bdes + eps), 'LineWidth', 2)
    % plot(ang, pow2db(Bmmse/max(Bmmse)), 'LineWidth', 2)
    % 
    % grid on
    % xlabel('Elevation (deg)')
    % ylabel('(dB)')
    % legend('Desired', 'CVX Covariance', 'Location', 'southoutside', 'Orientation', 'horizontal')
    % ylim([-20 1])
    % % xlim([0 90])
    % title('Transmit Beam Pattern')
end