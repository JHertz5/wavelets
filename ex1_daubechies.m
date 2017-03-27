%% Setup
clc
close all
clear variables

resolution = 64; maxTime = 32;
signalLength = 64*32;

time = (1:signalLength)./resolution;

%% Compute and plot Daubechies scaling function

% want to reproduce polynomials with degree 3
% scaling function must produce wavelets with 3+1=4 vanishing moments
% dbN where N is the number of vanishing moments
% db4 is selected as it has 4 vanishing moments

phi = zeros(1,signalLength);
[phi_T,~,~] = wavefun('db4',6); 
phi(1:length(phi_T))=phi_T;
figure
plot(time(1:length(phi_T)),phi_T,'LineWidth',2)
axis tight
xlabel('time \s')
ylabel('db4 scaling function')

%% Set Up Coefficient Computation

% initially found support to be 7, supported by .. saying it was 2N-1

m_degree = 0:4;
n_numCoefficients = 32;

t = ones(4,signalLength); % t^0 and initialise rest of matrix
t(2,:) = (0:signalLength-1)./64; % t^1
t(3,:) = t(2,:).^2; % t^2
t(4,:) = t(2,:).^3; % t^3
t(5,:) = t(2,:).^4; % t^4
% Construct the vectors of t^m, m = 0:degree-1 . 
% Also generate an extra t^m+1 to show failure
t = ones(m_degree(end),signalLength);
for order = m_degree+1
    t(order,:) = ((0:signalLength-1)./resolution).^(order - 1);
end

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

%% Reproduce the polynomials



t_reproduced = zeros(length(m_degree),signalLength);
error = t_reproduced;
for mIndex = m_degree+1
    title1 = sprintf('Reconstructed monomial, order = %i',mIndex-1);
    title2 = sprintf('Error, order = %i',mIndex-1);
    
    figure
    subplot(1,2,1)
    hold on
    
    % compute and plot reconstruction
    for nIndex = 1:n_numCoefficients
        shift = (nIndex-1) * resolution; % find the shift
        phiShifted = zeros(1,signalLength); % initialise phiShifted
        phiShifted((1 + shift):(length(phi_T) + shift)) = phi_T; % compute phiShifted
        phiShifted = phiShifted(1:signalLength); % crop to signalLength
        
        reproduction_contribution = c_coefficients(mIndex,nIndex) * phiShifted; % compute contribution at this m and n
        t_reproduced(mIndex,:) = t_reproduced(mIndex,:) + reproduction_contribution; % add to sum
        plot(time,reproduction_contribution)
    end
    h1 = plot(time,t_reproduced(mIndex,:),'b','LineWidth',2);
    h2 = plot(time,t(mIndex,:),'r--','LineWidth',2);
    legend([h1 h2],{'Reconstruction', 'Actual'}, 'Location', 'northwest')
    xlim([0 signalLength])
    xlabel('time \s')
    axis tight
    title(title1)
    
    % compute and plot error
    error(mIndex,:) = t(mIndex,:) - t_reproduced(mIndex,:); 
    
    subplot(1,2,2)
    plot(time,error(mIndex,:))
    xlabel('time \s')
    axis([0 32 -1 1]) 
    title(title2)
end
 
save('reproductionCoefficients.mat', 'c_coefficients')