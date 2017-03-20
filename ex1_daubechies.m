%% Setup
clc
close all

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
figure(1)
plot(phi,'r')
xlim([0,signalLength])
hold on
line(xlim,[0,0],'Linestyle','--')

%% Compute coefficients

firstZeroIndex = find(phi(2:end) == 0,1);
support = ceil(firstZeroIndex/resolution);

m_degree = 0:3;
n_numCoefficients = 0:32-support;

t = ones(4,signalLength);
%t(1,:) = ones(1, signalLength);
t(2,:) = (0:signalLength-1)./64;
t(3,:) = t_1.^2;
t(4,:) = t_1.^3;

figure(2)
hold on
xlim([0,signalLength/10]) % /by 10 to see plots clearly
plot(t(1,:))
plot(t(2,:))
plot(t(3,:))
plot(t(4,:))
hold off

c_coefficients = zeros(length(m_degree),length(n_numCoefficients));
for mIndex = m_degree+1
    for nIndex = n_numCoefficients+1
        shift = (nIndex-1) * resolution;
        phiShifted = zeros(1,signalLength);
        phiShifted((1 + shift):(support*resolution + shift)) = phi(1:support*resolution);
        
        c_coefficients(mIndex, nIndex) = dot(t(mIndex,:),phiShifted);
    end
end

t_reproduced = zeros(length(m_degree),signalLength);
for mIndex = m_degree+1
    figure(mIndex+3)
    hold on
    for nIndex = n_numCoefficients+1
        shift = (nIndex-1) * resolution;
        phiShifted = zeros(1,signalLength);
        phiShifted((1 + shift):(support*resolution + shift)) = phi(1:support*resolution);
        
        reproduction_contribution = c_coefficients(mIndex,nIndex) * phiShifted;
        t_reproduced(mIndex,:) = t_reproduced(mIndex,:) + reproduction_contribution;
        plot(reproduction_contribution,':')
    end
    plot(t_reproduced(mIndex,:),'b','LineWidth',2)
end
 
