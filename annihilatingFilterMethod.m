function [ h_annihilatingFilter, tk_locations, ak_weights, y  ] = annihilatingFilterMethod(tau_moments)
%annihilatingFilterMethod
% calcualte annihilating filter values for signal
% and return locations and weights

K = 2;
N = length(tau_moments)-1;

%% Find annihilating filter by solving equation

h_annihilatingFilter = zeros(K+1,1);
h_annihilatingFilter(1) = 1;

eqn_tauVect1 = tau_moments(K+1:N+1); % column vector of tau from K to N
eqn_tauMatrix = zeros(N-K+1, K); % initialise N-K x K-1 matrix
for index = 0:K-1
    eqn_tauMatrix(:,K-index) = tau_moments(1+index : N-K+1+index);
end

% solve eqn_tauMatrix * eqn_hVect = eqn_tauVect1;
%eqn_hVect = eqn_tauMatrix^-1 * -eqn_tauVect1;
eqn_hVect = eqn_tauMatrix \ -eqn_tauVect1;

h_annihilatingFilter(2:end) = eqn_hVect; % solve equation to find filter values

%% Find convolution

y = conv(tau_moments,h_annihilatingFilter);

%% Find locations

tk_locations = roots(h_annihilatingFilter); % locations are roots of filter

% Print location values
%disp('Locations:')
%disp(tk_locations.')

%% Solve Vandermonde system to find weights

eqn_locationsMatrix = [ 1 1; tk_locations(1) tk_locations(2) ];
eqn_tauVect2 = tau_moments(1:K);

ak_weights = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

end

