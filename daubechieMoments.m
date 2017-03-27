function [ moments, phi, coefficients ] = daubechieMoments( xt, degree )
%DAUBECHIEMOMENTS Creates moments of function with Daubechie of degree.
%   Makes Q7 easier with different N
    
    % Generate the Daubechie scaling function
    signalLength = length(xt);
    waveorder = num2str(degree);
    wavetype = ['db' waveorder];
    phi = zeros(1, signalLength);
    [phiT, ~, ~] = wavefun(wavetype, 6);
    phi(1:length(phiT))=phiT;
    
    % Constants for later
    L = 2 * degree - 1;
    nVec = 0:32-L;
    T = 64;
    
    % Create t
    % Construct the vectors of t^m, m = 0:degree-1 . 
    tVals = ones(signalLength, degree);
    for order = 2:degree
        tVals(:, order) = ((0:signalLength-1)/T).^(order - 1);
    end
    
    % Create shifted phi
    allPhi = zeros(length(nVec), length(phi));
    for n = nVec
        allPhi(n+1, :) = [zeros(1, n*T) phi(1:end - n*T)];
    end
    
    % Acquire coefficients
    coefficients =  (1/T) * allPhi * tVals;
    
    % Take xt and make yn
    yn = allPhi * xt;
    % Moments
    moments = yn' * coefficients;

end

