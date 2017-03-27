%% Setup
clc
close all
clear variables

K = 2;

resolution = 64; maxTime = 32;
signalLength = 64*32;

%% Create stream of Diracs

ak_weights = [7; 4];
tk_locations = [7.5; 23];

% initialise vector and add diracs
x_diracsStream = zeros(1,signalLength); 
x_diracsStream(tk_locations(1)*resolution+1) = ak_weights(1);
x_diracsStream(tk_locations(2)*resolution+1) = ak_weights(2);

figure
stem(x_diracsStream,'x')
xlim([0 signalLength])

disp('Original Dirac values -')
disp('Locations:')
disp(tk_locations')
disp('Weights:')
disp(ak_weights')

%% Daubechies
N_range = 5:8;
noiseSigmas = sqrt([0.001 0.01 0.1 1 10]);
% Multiplying by a gives variance of a^2. 

for order = N_range
    % Create moments
    [moments, phi, ~] = momentsGenerator(x_diracsStream, order);
    % Duplicate
    momentsNoise = repmat(moments.', length(noiseSigmas), 1);
    
    % Add different values of noise.
    for index = 1:length(noiseSigmas)
        momentsNoise(index,:) = noiseSigmas(index) * randn(1, length(moments)) + momentsNoise(index,:);
    end
    
    % Picture time!
    figure('position',[0 0 1280 800]);

    subplot(3, 2, 1);
    plot(phi);
    axis([0 7*resolution -0.4 1.2]);
    the_title = ['Order ' num2str(order) ' Daubechies Sampling Kernel'];
    title(the_title);
    xlabel('Time');
    
    for index = 1:length(noiseSigmas)
        subplot(3, 2, index+1);
        stem(momentsNoise(index, :),'x');
        xlabel('m');
        title(['Moment, \sigma^2: ' num2str(noiseSigmas(index)^2)]);
    end
    
    filename = ['noisyMoments' num2str(order)];
    save(filename, 'momentsNoise', 'phi');
end
