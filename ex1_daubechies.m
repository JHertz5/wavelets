%% Setup
clc
close all
clear variables

resolution = 64; maxTime = 32;
signalLength = 64*32;

% want to reproduce polynomials with degree 3
% scaling function must produce wavelets with 3+1=4 vanishing moments
% dbN where N is the number of vanishing moments
% db4 is selected as it has 4 vanishing moments

%% Compute and plot Daubechies scaling function

phi = zeros(1,signalLength);
[phi_T, psi_T, xval] = wavefun('db4',6); 
phi(1:length(phi_T))=phi_T;
figure
plot(phi_T,'LineWidth',2)
xlim([0 length(phi_T)])
xlabel('time (s/64)')
ylabel('db4 scaling function')

%% Compute coefficients

firstZeroIndex = find(phi(2:end) == 0,1); % find first zero after position 1
support = ceil(firstZeroIndex/resolution); % compute support

m_degree = 0:3;
n_numCoefficients = 0:32-support;

t = ones(4,signalLength); % t^0 and initialise rest of matrix
t(2,:) = (0:signalLength-1)./64; % t^1
t(3,:) = t(2,:).^2; % t^2
t(4,:) = t(2,:).^3; % t^3

%% Compute coefficients

c_coefficients = zeros(length(m_degree),length(n_numCoefficients));
for mIndex = m_degree+1
    for nIndex = n_numCoefficients+1
        shift = (nIndex-1) * resolution; % find the shift
        phiShifted = zeros(1,signalLength); % initialise phiShifted
        phiShifted((1 + shift):(support*resolution + shift)) = phi(1:support*resolution); % compute phiShifted
                
        c_coefficients(mIndex, nIndex) = dot(t(mIndex,:),phiShifted)./resolution; % inner product to find c
    end
end

%% Reproduce the polynomials

t_reproduced = zeros(length(m_degree),signalLength);
for mIndex = m_degree+1
    figure
    hold on
    for nIndex = n_numCoefficients+1
        shift = (nIndex-1) * resolution; % find the shift
        phiShifted = zeros(1,signalLength); % initialise phiShifted
        phiShifted((1 + shift):(support*resolution + shift)) = phi(1:support*resolution); % compute phiShifted
        
        reproduction_contribution = c_coefficients(mIndex,nIndex) * phiShifted; % compute contribution at this m and n
        t_reproduced(mIndex,:) = t_reproduced(mIndex,:) + reproduction_contribution; % add to sum
        plot(reproduction_contribution)
    end
    h1 = plot(t_reproduced(mIndex,:),'b','LineWidth',2);
    h2 = plot(t(mIndex,:),'r--','LineWidth',2);
    legend([h1 h2],{'Reconstruction', 'Actual'}, 'Location', 'northwest')
    xlim([0 signalLength])
    xlabel('time (s/64)')
end
 
save('reproductionCoefficients.mat', 'c_coefficients')