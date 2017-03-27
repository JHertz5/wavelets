%% Setup
clc
close all
clear variables

resolution = 64; maxTime = 32;
signalLength = 64*32;

time = (1:signalLength)./resolution;

%% Find the Dual basis of the B-Spline 
syms z;
%construct P
subPoly = 2 - z - z^-1;
sum_P = 1 + subPoly + (5/8)*subPoly^2 + (5/16)*subPoly^3;
% comes from H(z) = P(z)/G(z)
H = expand(((z+1)^4) * (sum_P));
%sum2poly doesn't work on negative indices, so shift H up by 3
dual_phi_T = sym2poly((z^3)*H);

%% Calculate the B spline scaling function
[beta0, psi_T, xval] = wavefun('haar', 1);
phi_T = conv(beta0, conv(beta0,conv(beta0,conv(beta0, beta0))));
phi_T = phi_T(4:end-2); % cut off zero parts of function

figure
plot(phi_T,'LineWidth',2)
axis tight
xlabel('Tap')
ylabel('B3 Scaling Function')

%{
subplot(1,2,2)
plot(dual_phi_T,'LineWidth',2)
axis tight
xlabel('Tap')
title('Dual Basis of B3 Scaling Function')

%% Set Up Coefficient Computation

m_degree = 0:3;
n_numCoefficients = 1024;

t = ones(4,signalLength); % t^0 and initialise rest of matrix
t(2,:) = (0:signalLength-1)./64; % t^1
t(3,:) = t(2,:).^2; % t^2
t(4,:) = t(2,:).^3; % t^3
t(5,:) = t(2,:).^4; % t^4

%% Compute coefficients

c_coefficients = zeros(length(m_degree),n_numCoefficients);
for mIndex = m_degree+1
    for nIndex = 1:n_numCoefficients
        shift = (nIndex-1) * 2;% resolution; % find the shift
        phiShifted = zeros(1,signalLength); % initialise phiShifted
        phiShifted((1 + shift):(length(phi_T) + shift)) = dual_phi_T; % compute phiShifted
        phiShifted = phiShifted(1:signalLength); % crop to signalLength
                
        c_coefficients(mIndex, nIndex) = dot(t(mIndex,:),phiShifted)./(resolution*4); % inner product to find c
    end
end

%% Reproduce the polynomials

t_reproduced = zeros(length(m_degree),signalLength);
figure

for mIndex = m_degree+1
    title1 = sprintf('Reconstructed monomial, order = %i',mIndex-1);    
    
    subplot(2,2,mIndex)
    hold on
    
    % compute and plot reconstruction
    for nIndex = 1:n_numCoefficients
        shift = (nIndex-1) * 2;%resolution; % find the shift
        phiShifted = zeros(1,signalLength); % initialise phiShifted
        phiShifted((1 + shift):(length(phi_T) + shift)) = phi_T; % compute phiShifted
        phiShifted = phiShifted(1:signalLength); % crop to signalLength
        
        reproduction_contribution = c_coefficients(mIndex,nIndex) * phiShifted; % compute contribution at this m and n
        t_reproduced(mIndex,:) = t_reproduced(mIndex,:) + reproduction_contribution; % add to sum
        %plot(time,reproduction_contribution)
        plot(reproduction_contribution)
    end
    %h1 = plot(time,t_reproduced(mIndex,:),'b','LineWidth',2);
    %h2 = plot(time,t(mIndex,:),'r--','LineWidth',2);
    h1 = plot(t_reproduced(mIndex,:),'b','LineWidth',2);
    h2 = plot(t(mIndex,:),'r--','LineWidth',2);
    legend([h1 h2],{'Reconstruction', 'Actual'}, 'Location', 'northwest')
    xlim([0 signalLength])
    xlabel('time \s')
    axis tight
    ylabel(title1)
end
%}