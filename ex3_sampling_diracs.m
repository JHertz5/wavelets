%% Setup
clc
close all
clear variables

K = 2;

resolution = 64; maxTime = 32;
signalLength = 64*32;

%% Create stream of Diracs

ak_weights = [1.24 0.8];
tk_locations = [9.25 20.4375];

x_diracsStream = zeros(1,signalLength);
x_diracsStream(tk_locations(1)*resolution) = ak_weights(1);
x_diracsStream(tk_locations(2)*resolution) = ak_weights(2);

figure
stem(x_diracsStream,'x')
xlim([0 signalLength])

disp('Original Dirac values -')
disp('Locations:')
disp(tk_locations.')
disp('Weights:')
disp(ak_weights.')

%% Compute and plot Daubechies scaling function

% N >= 2K-1 (p70)
% N = 2*(2)-1 = 3 
% N+1 = 3+1 = 4 -> db4

phi = zeros(1,signalLength);
[phi_T, psi_T, xval] = wavefun('db4',6); 
phi(1:length(phi_T))=phi_T;

%% Sample signal using Daubechies scaling function

firstZeroIndex = find(phi(2:end) == 0,1); % find first zero after position 1
support = ceil(firstZeroIndex/resolution); % compute support

m_degree = 0:3;
n_numSamples = 32-support;

y_samples = zeros(1, n_numSamples);
for sampleIndex = 1:n_numSamples
    shift = sampleIndex * resolution; % find the shift
    phiShifted = zeros(1,signalLength); % initialise phiShifted
    phiShifted((1 + shift):(support*resolution + shift)) = phi(1:support*resolution); % compute phiShifted
    
    y_samples(sampleIndex) = dot(x_diracsStream,phiShifted);
end
figure
stem(y_samples)

