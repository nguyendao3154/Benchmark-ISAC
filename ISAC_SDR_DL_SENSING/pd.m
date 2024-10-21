% MATLAB code to plot SNR vs Probability of Detection using LRT

% clear; clc; close all;

% Parameters
matrix = ones(1,201);
% for a = 3:6
 a = 4;
P_FA = 1*10^-a;  % Probability of False Alarm (fixed)
SNR_dB = 0:0.1:20;  % SNR values in dB
SNR_linear = 10.^(SNR_dB/10);  % Convert SNR to linear scale

% Calculate the threshold gamma based on P_FA (for Gaussian noise)
gamma = sqrt(2) * erfcinv(2 * P_FA);  % erfcinv is the inverse complementary error function
test = sqrt(-log(P_FA));
% Calculate Probability of Detection (P_D) for each SNR value
P_D = zeros(size(SNR_dB));  % Initialize P_D vector

for i = 1:length(SNR_dB)
    % Q-function for Gaussian noise and signal
    P_D(i) = 0.5 * erfc((gamma - sqrt(SNR_linear(i))) / sqrt(2));
end

% Plot SNR vs Probability of Detection
    hold on;
    plot(SNR_dB, P_D, 'LineWidth', 2);
    grid on;
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('Probability of Detection, P_D', 'FontSize', 12);
    title('SNR vs Probability of Detection (LRT)', 'FontSize', 14);
    ylim([0 1]); 
    
    % Display P_FA on plot
    legend(['P_{FA} = ' num2str(P_FA)], 'Location', 'southeast'); 
    matrix = [matrix; P_D]
% end
matrix(1,:) = 0:0.1:20;
% writematrix(transpose(matrix), 'P_D.csv');