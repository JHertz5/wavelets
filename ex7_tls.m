%% Setup
clear;
close all;
clc;

% Some constants
resolution = 64;
maxTime = 32;
signalLength = resolution * maxTime;
numberOfIterations = 6;
K = 2;
time = (1:signalLength)./resolution;

%% Create stream of Diracs

ak_weights = [7; 4];
tk_locations = [7.5; 23];

% initialise vector and add diracs
x_diracsStream = zeros(1,signalLength); 
x_diracsStream(tk_locations(1)*resolution+1) = ak_weights(1);
x_diracsStream(tk_locations(2)*resolution+1) = ak_weights(2);

%% Acquire moments and Calculate Diracs
for moment = 5:8
    if moment == 5
        load('noisyMoments5.mat');
    elseif moment == 6
        load('noisyMoments6.mat');
    elseif moment == 7
        load('noisyMoments7.mat');
    elseif moment == 8
        load('noisyMoments8.mat');
    end
    noiseSigmas = size(momentsNoise, 1);

    % Calculate the Diracs
    for noiseIndex = 1:noiseSigmas
        tau_moments = momentsNoise(noiseIndex, :);
        [h, tk_locations_est, ak_weights_est] = annihilatingFilterMethodTLS(tau_moments);

            x_diracsStream_est = zeros(1,signalLength);
            for diracIndex = 1:K
                x_diracsStream_est(uint32(tk_locations_est(diracIndex)*resolution+1)) = ak_weights_est(diracIndex);
            end

            figure

            stem(time, x_diracsStream_est,'x');
            axis tight
            title('Reconstructed Signal');
            xlabel('time');
    end
end
