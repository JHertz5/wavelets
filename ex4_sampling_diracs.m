%% Setup
clc
close all
clear variables

load reproductionCoefficients.mat

K = 2;

resolution = 64; maxTime = 32;
signalLength = 64*32;

%% Create stream of Diracs

ak_weights = [5; 2];
tk_locations = [12.5; 13];

% initialise vector and add 
x_diracsStream = zeros(1,signalLength); 
x_diracsStream(tk_locations(1)*resolution) = ak_weights(1);
x_diracsStream(tk_locations(2)*resolution) = ak_weights(2);

figure
stem(x_diracsStream,'x')
xlim([0 signalLength])

disp('Original Dirac values -')
disp('Locations:')
disp(tk_locations)
disp('Weights:')
disp(ak_weights)

%% Compute and plot Daubechies scaling function

% N >= 2K-1 (p70)
% N = 2*(2)-1 = 3 
% N+1 = 3+1 = 4 -> db4

phi = zeros(1,signalLength);
[phi_T,~,~] = wavefun('db4',6); 
phi(1:length(phi_T))=phi_T;

%% Sample signal using Daubechies scaling function

firstZeroIndex = find(phi(2:end) == 0,1); % find first zero after position 1
support = 7; %ceil(firstZeroIndex/resolution); % compute support

m_degree = 0:3;
n_numSamples = 32-support;

y_samples = zeros(1, n_numSamples);
for sampleIndex = 1:n_numSamples+1
    shift = (sampleIndex-1) * resolution; % find the shift
    phiShifted = zeros(1,signalLength); % initialise phiShifted
    phiShifted((1 + shift):(length(phi_T) + shift)) = phi_T; % compute phiShifted
    phiShifted = phiShifted(1:signalLength); % 'crop' phiShifted

    y_samples(sampleIndex) = dot(x_diracsStream,phiShifted);
end

% plot samples
figure
stem(y_samples)

%% Retrieve N+1 moments of signal

s_moments = zeros(length(m_degree),1);
for degreeIndex = m_degree+1
    s_moments(degreeIndex) = dot(c_coefficients(degreeIndex,:),y_samples);
end

%% Apply annihilating filter method

[tk_locations_est, ak_weights_est] = annihilatingFilterMethod(s_moments,true);

tk_locations_est = tk_locations_est + 0.0156; % hack

disp('Estimated Dirac values -')
disp('Locations:')
disp(tk_locations_est)
disp('Weights:')
disp(ak_weights_est)