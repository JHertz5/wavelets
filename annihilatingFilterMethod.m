function [ tk_locations, ak_weights  ] = annihilatingFilterMethod( tau_signal, showPlots )
%annihilatingFilterMethod
% calcualte annihilating filter values for signal
% and return locations and weights

K = 2;
N = size(tau_signal,2)-1;

%% Find annihilating filter by solving equation

h_annihilatingFilter = zeros(K+1,1);
h_annihilatingFilter(1) = 1;

eqn_tauVect1 = tau_signal(K+1:N+1).'; % column vector of tau from K to N
eqn_tauMatrix = zeros(N-K+1, K); % initialise N-K x K-1 matrix
for index = 0:K-1
    eqn_tauMatrix(:,K-index) = tau_signal(1+index : N-K+1+index);
end

% solve eqn_tauMatrix * eqn_hVect = eqn_tauVect1;
%eqn_hVect = eqn_tauMatrix^-1 * -eqn_tauVect1;
eqn_hVect = eqn_tauMatrix \ -eqn_tauVect1;

h_annihilatingFilter(2:end) = eqn_hVect; % solve equation to find filter values

%% Plot convolution

y = conv(tau_signal,h_annihilatingFilter);
if showPlots
    % Plot tau
    subplot(311)
    stem(tau_signal,'x')
    ylabel('Tau')

    % Plot filter
    subplot(312)
    stem(h_annihilatingFilter,'x')
    ylabel('Annihilating Filter')

    % Plot convolution - only middle two values need to be 0
    subplot(313)
    stem(y,'x');
    ylabel('Convolution')
end

%% Find locations

tk_locations = roots(h_annihilatingFilter); % locations are roots of filter

% Print location values
disp('Locations:')
disp(tk_locations.')

%% Solve Vandermonde system to find weights

eqn_locationsMatrix = [ 1 1; tk_locations(1) tk_locations(2) ];
eqn_tauVect2 = tau_signal(1:K).';

ak_weights = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

% Print weight values
disp('Weights:')
disp(ak_weights.')

end

