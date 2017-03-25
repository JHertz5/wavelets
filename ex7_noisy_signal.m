%% Setup
clc
close all
clear variables

load noisyMoments.mat % sets K and N as well


disp('Correct Values')
% Print location values
disp('Locations:')
disp(tk_locations.')
% Print weight values
disp('Weights:')
disp(ak_weights.')

%% Perform Total Least-Squares approach to find annihilating filter

% Fill in the matrix row by row
toeplitz = zeros(K+1,N-K);
for index = 1:K+1
    toeplitz(index,:) = flip(s_noisyMoments(index:K+index));
end

[U, lambda, V] = svd(toeplitz);
 
h_annihilatingFilter = V(:,end); % solve equation to find filter values

%% Plot convolution

y = conv(s_noisyMoments,h_annihilatingFilter);

% Plot tau
subplot(311)
stem(s_noisyMoments,'x')
ylabel('Noisy Moments')

% Plot filter
subplot(312)
stem(h_annihilatingFilter,'x')
ylabel('Annihilating Filter')

% Plot convolution - only middle two values need to be 0
subplot(313)
stem(y,'x');
ylabel('Convolution')

%% Find locations

tk_locations_est = roots(h_annihilatingFilter); % locations are roots of filter

%% Solve Vandermonde system to find weights

eqn_locationsMatrix = [ 1 1; tk_locations_est(1) tk_locations_est(2) ];
eqn_tauVect2 = s_noisyMoments(1:K);

ak_weights_est = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

%% Print LTS results

disp('LTS only')
% Print location values
disp('Locations:')
disp(tk_locations_est.')
% Print weight values
disp('Weights:')
disp(ak_weights_est.')





%% Perform Cadzow routine
% iteration 1, using lambda from earlier
diagLambda = diag(lambda);

for index = 1:length(diagLambda)-K
    [noiseCoefficient, noiseIndex] = min(diagLambda);
    diagLambda(noiseIndex) = 0;
end

lambda_denoised1 = diag(diagLambda);
toeplitz_denoised1 = U*lambda_denoised1*V';
tau_momentsDenoised1 = [toeplitz_denoised1(1,3); toeplitz_denoised1(2,3); toeplitz_denoised1(3,3); toeplitz_denoised1(3,2); toeplitz_denoised1(3,1)];

%% Perform denoised Total Least-Squares approach to find annihilating filter

[~, ~, V_denoised1] = svd(toeplitz_denoised1);
 
h_annihilatingFilter_denoised1 = V_denoised1(:,end); % solve equation to find filter values

%% Find locations

tk_locations_estCadzow1 = roots(h_annihilatingFilter_denoised1); % locations are roots of filter

%% Solve Vandermonde system to find weights

eqn_locationsMatrix = [ 1 1; tk_locations_estCadzow1(1) tk_locations_estCadzow1(2) ];
eqn_tauVect2 = tau_momentsDenoised1(1:K);

ak_weights_estCadzow1 = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

%% Print LTS results

disp('Cadzow + LTS')
% Print location values
disp('Locations:')
disp(tk_locations_estCadzow1.')
% Print weight values
disp('Weights:')
disp(ak_weights_estCadzow1.')







%% iteration 2+

[U_denoised1, lambda_denoised1, V_denoised1] = svd(toeplitz_denoised1);
diagLambda = diag(lambda_denoised1);

for index = 1:length(diagLambda)-K
    [noiseCoefficient, noiseIndex] = min(diagLambda);
    diagLambda(noiseIndex) = 0;
end

lambda_denoised2 = diag(diagLambda);
toeplitz_denoised2 = U_denoised1*lambda_denoised2*V_denoised1';

tau_momentsDenoised2 = [toeplitz_denoised2(1,3); toeplitz_denoised2(2,3); toeplitz_denoised2(3,3); toeplitz_denoised2(3,2); toeplitz_denoised2(3,1)];

%% Perform denoised Total Least-Squares approach to find annihilating filter

[~, ~, V_denoised2] = svd(toeplitz_denoised1);
 
h_annihilatingFilter_denoised2 = V_denoised2(:,end); % solve equation to find filter values

%% Find locations

tk_locations_estCadzow2 = roots(h_annihilatingFilter_denoised2); % locations are roots of filter

%% Solve Vandermonde system to find weights

eqn_locationsMatrix = [ 1 1; tk_locations_estCadzow2(1) tk_locations_estCadzow2(2) ];
eqn_tauVect2 = tau_momentsDenoised2(1:K);

ak_weights_estCadzow2 = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

%% Print LTS results

disp('Cadzow (x2) + LTS')
% Print location values
disp('Locations:')
disp(tk_locations_estCadzow2.')
% Print weight values
disp('Weights:')
disp(ak_weights_estCadzow2.')

