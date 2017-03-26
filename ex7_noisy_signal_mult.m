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

%% Use annihilating filter to find locations and weights

tk_locations_est = roots(h_annihilatingFilter); % locations are roots of filter

eqn_locationsMatrix = [ 1 1; tk_locations_est(1) tk_locations_est(2) ];
eqn_tauVect2 = s_noisyMoments(1:K);

ak_weights_est = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

disp('LTS only')
% Print location values
disp('Locations:')
disp(tk_locations_est.')
% Print weight values
disp('Weights:')
disp(ak_weights_est.')


lambda_de = lambda;
U_de = U;
V_de = V;
%% Perform Cadzow routine
% iteration 1, using lambda from earlier
iterations = 10;
for iterationIndex = 1:iterations

    diagLambda = diag(lambda_de);

    for index = 1:length(diagLambda)-K
        [noiseCoefficient, noiseIndex] = min(diagLambda);
        diagLambda(noiseIndex) = 0;
    end

    lambda_de = diag(diagLambda);
    toeplitz_de = U_de*lambda_de*V_de';

    if(size(toeplitz_de) == [3 3])
        %first diag
        ave = (toeplitz_de(1,2) + toeplitz_de(2,3))/2;
        toeplitz_de(1,2) = ave;
        toeplitz_de(2,3) = ave;

        %second diag
        ave = (toeplitz_de(1,1) + toeplitz_de(2,2) + toeplitz_de(3,3))/3;
        toeplitz_de(1,1) = ave;
        toeplitz_de(2,2) = ave;
        toeplitz_de(3,3) = ave;

        %third diag
        ave = (toeplitz_de(2,1) + toeplitz_de(3,2))/2;
        toeplitz_de(2,1) = ave;
        toeplitz_de(3,2) = ave;
    else
        disp('TOEPLITZ DIAG AVERAGER DOESNT WORK FOR THIS SIZE TOEPLITZ')
    end

end
%% Perform denoised Total Least-Squares approach to find annihilating filter

[~, ~, V_de1] = svd(toeplitz_de);
h_annihilatingFilter_de1 = V_de1(:,end); % solve equation to find filter values

%% Use annihilating filter to find locations and weights

s_momentsDenoised1 = [toeplitz_de(1,3); toeplitz_de(2,3); toeplitz_de(3,3); toeplitz_de(3,2); toeplitz_de(3,1)];

tk_locations_estCadzow1 = roots(h_annihilatingFilter_de1); % locations are roots of filter

eqn_locationsMatrix = [ 1 1; tk_locations_estCadzow1(1) tk_locations_estCadzow1(2) ];
eqn_tauVect2 = s_momentsDenoised1(1:K);

ak_weights_estCadzow1 = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

disp('Cadzow + LTS')
% Print location values
disp('Locations:')
disp(tk_locations_estCadzow1.')
% Print weight values
disp('Weights:')
disp(ak_weights_estCadzow1.')
