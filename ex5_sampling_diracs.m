%% Setup
clc
close all
clear variables

load reproductionCoefficients.mat
load project_files_&_data/samples.mat

K = 2;

resolution = 64; maxTime = 32;
signalLength = 64*32;
time = (1:signalLength)./resolution;

%% Retrieve N+1 moments of signal

m_degree = 0:3;

s_moments = zeros(length(m_degree),1);
for degreeIndex = m_degree+1
    s_moments(degreeIndex) = dot(c_coefficients(degreeIndex,:),y_sampled);
end

%% Apply annihilating filter method

[h, tk_locations_est, ak_weights_est, y] = annihilatingFilterMethod(s_moments);

%% Display and Plot Results

% initialise vector and add diracs
x_diracsStream_est = zeros(1,signalLength); 
x_diracsStream_est(uint32(tk_locations_est(1)*resolution+1)) = ak_weights_est(1);
x_diracsStream_est(uint32(tk_locations_est(2)*resolution+1)) = ak_weights_est(2);

stem(time,x_diracsStream_est,'x')
axis tight
ylabel('amplitude')
xlabel('time /s')
title('Reconstruction of Dirac Stream')

disp('Estimated Dirac values -')
disp('Locations:')
disp(tk_locations_est)
disp('Weights:')
disp(ak_weights_est)

%% Compute and plot Daubechies scaling function

% N >= 2K-1 (p70)
% N = 2*(2)-1 = 3 
% N+1 = 3+1 = 4 -> db4

phi = zeros(1,signalLength);
[phi_T,~,~] = wavefun('db4',6); 
phi(1:length(phi_T))=phi_T;

%% Sample signal using Daubechies scaling function

m_degree = 0:3;
n_numSamples = 32;

y_sampled2 = zeros(1, n_numSamples);
for sampleIndex = 1:n_numSamples
    shift = (sampleIndex-1) * resolution; % find the shift
    phiShifted = zeros(1,signalLength); % initialise phiShifted
    phiShifted((1 + shift):(length(phi_T) + shift)) = phi_T; % compute phiShifted
    phiShifted = phiShifted(1:signalLength); % crop to signalLength

    y_sampled2(sampleIndex) = dot(x_diracsStream_est,phiShifted); % compute samples
end

figure

subplot(2,1,1)
stem(y_sampled)
axis tight
ylabel('amplitude')
xlabel('time /s')
title('Original Samples')

subplot(2,1,2)
stem(y_sampled2)
axis tight
ylabel('amplitude')
xlabel('time /s')
title('Samples of Reconstruction')