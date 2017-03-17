clc

load project_files_&_data/tau.mat

K = 2;
N = size(tau,2)-1;

h_annihilatingFilter = zeros(K+1,1);
h_annihilatingFilter(1) = 1;

%% Find annihilating filter by solving equation

eqn_tauVect = tau(K+1:N+1).'; % column vector of tau from K to N
eqn_tauMatrix = zeros(N-K+1, K);
for index = 0:K-1
    eqn_tauMatrix(:,K-index) = tau(1+index : N-K+1+index);
end

% solve eqn_tauMatrix * eqn_hVect = eqn_tauVect;
%eqn_hVect = eqn_tauMatrix^-1 * -eqn_tauVect;
eqn_hVect = eqn_tauMatrix \ -eqn_tauVect;

h_annihilatingFilter(2:end) = eqn_hVect;

%% Plot convolution

y = conv(tau,h_annihilatingFilter);

% Plot tau
subplot(311)
ylabel('Tau')
stem(tau,'x')

% Plot filter
subplot(312)
ylabel('Annihilating Filter')
stem(h_annihilatingFilter,'x')

% Plot convolution
subplot(313)
ylabel('Convolution')
stem(y,'x')

%% Find locations and weights

tk_locations = roots(h_annihilatingFilter);

eqn_locationsMatrix = [ 1 1; tk_locations(1) tk_locations(2) ];
eqn_tauVect = tau(1:K).';

ak_weights = eqn_locationsMatrix \ eqn_tauVect;