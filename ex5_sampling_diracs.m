%% Setup
clc
close all
clear variables

load reproductionCoefficients.mat
load project_files_&_data/samples.mat

K = 2;

resolution = 64; maxTime = 32;
signalLength = 64*32;

%% Retrieve N+1 moments of signal

m_degree = 0:3;

s_moments = zeros(length(m_degree),1);
for degreeIndex = m_degree+1
    s_moments(degreeIndex) = dot(c_coefficients(degreeIndex,:),y_sampled);
end

%% Apply annihilating filter method

[tk_locations_est, ak_weights_est] = annihilatingFilterMethod(s_moments,true);

disp('Estimated Dirac values -')
disp('Locations:')
disp(tk_locations_est)
disp('Weights:')
disp(ak_weights_est)