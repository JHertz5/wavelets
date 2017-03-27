%% Setup
clc
close all
clear variables

load project_files_&_data/tau.mat

tau = tau';

resolution = 64; maxTime = 32;
signalLength = 64*32;
time = (1:signalLength)./resolution;

%% Apply Annihilating Filter Method

[ h, tk_locations_est, ak_weights_est, y ] = annihilatingFilterMethod(tau);

% initialise vector and add diracs
x_diracsStream_est = zeros(1,signalLength); 
x_diracsStream_est(uint32(tk_locations_est(1)*resolution+1)) = ak_weights_est(1);
x_diracsStream_est(uint32(tk_locations_est(2)*resolution+1)) = ak_weights_est(2);

figure
stem(time,x_diracsStream_est,'x')
axis tight
ylabel('amplitude')
xlabel('time /s')
%title('Reconstruction of Dirac Stream')

disp('Estimated Dirac values -')
disp('Locations:')
disp(tk_locations_est)
disp('Weights:')
disp(ak_weights_est)