clc

% want to reproduce polynomials with degree 3
% scaling function must produce wavelets with 3+1=4 vanishing moments
% dbN where N is the number of vanishing moments
% db4 is selected as it has 4 vanishing moments

resolution = 64; maxTime = 32;
signalLength = 64*32;

phi = zeros(1,signalLength);
[phi_T, psi_T, xval] = wavefun('db4',6); 
phi(1:length(phi_T))=phi_T;
plot(phi,'r')
xlim([0,signalLength])
hold on
line(xlim,[0,0],'Linestyle','--')

firstZeroIndex = find(phi(2:end) == 0,1);
support = ceil(firstZeroIndex/resolution);

m_degree = 0:3;
n_numCoefficients = 0:32-support;

phi_dual = inv(phi).';
