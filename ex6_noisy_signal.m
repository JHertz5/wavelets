%% Setup
clc
close all
clear variables

K = 2;

resolution = 64; maxTime = 32;
signalLength = 64*32;

%% Create stream of Diracs

ak_weights = [5; 2];
tk_locations = [12; 23];

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

%% Compute and plot Daubechies scaling function

% N > 2K -> N > 2K
% N = 2*(2)-1 = 3 
% N+1 = 3+1 = 4 -> db4

N = 5;

phi = zeros(1,signalLength);
[phi_T,~,~] = wavefun('db6',6); 
phi(1:length(phi_T))=phi_T;
plot(phi)

%% Set Up Coefficient Computation

m_degree = 0:N-1;
n_numCoefficients = 32;

t = ones(5,signalLength); % t^0 and initialise rest of matrix
t(2,:) = (0:signalLength-1)./64; % t^1
t(3,:) = t(2,:).^2; % t^2
t(4,:) = t(2,:).^3; % t^3
t(5,:) = t(2,:).^4; % t^4

%% Compute coefficients

c_coefficients = zeros(length(m_degree),n_numCoefficients);
for mIndex = m_degree+1
    for nIndex = 1:n_numCoefficients
        shift = (nIndex-1) * resolution; % find the shift
        phiShifted = zeros(1,signalLength); % initialise phiShifted
        phiShifted((1 + shift):(length(phi_T) + shift)) = phi_T; % compute phiShifted
        phiShifted = phiShifted(1:signalLength); % crop to signalLength
                
        c_coefficients(mIndex, nIndex) = dot(t(mIndex,:),phiShifted)./resolution; % inner product to find c
    end
end

%% Sample signal using Daubechies scaling function

n_numSamples = 32;

y_sampled = zeros(1, n_numSamples);
for sampleIndex = 1:n_numSamples
    shift = (sampleIndex-1) * resolution; % find the shift
    phiShifted = zeros(1,signalLength); % initialise phiShifted
    phiShifted((1 + shift):(length(phi_T) + shift)) = phi_T; % compute phiShifted
    phiShifted = phiShifted(1:signalLength); % crop to signalLength

    y_sampled(sampleIndex) = dot(x_diracsStream,phiShifted); % compute samples
end

% plot samples
figure
stem(y_sampled)

%% Retrieve N+1 moments of signal

s_moments = zeros(length(m_degree),1);
for degreeIndex = m_degree+1
    s_moments(degreeIndex) = dot(c_coefficients(degreeIndex,:),y_sampled);
end

%% Add noise to s_moments

sigmas = s_moments./1000; % setting standard deviation to be 10% of s_moments values
epsilon_noiseComponent = sigmas.*randn(size(s_moments,1), size(s_moments,2));
s_noisyMoments = s_moments + epsilon_noiseComponent;

save('noisyMoments.mat', 's_noisyMoments', 'K', 'N', 'ak_weights', 'tk_locations')
