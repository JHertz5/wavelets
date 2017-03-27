function [ moments, phi, coefficients ] = momentsGenerator( x_diracsStream, degree )
%DAUBECHIEMOMENTS Creates moments of function with Daubechie of degree.
%   Makes Q7 easier with different N
    
    % Generate the Daubechie scaling function
    signalLength = length(x_diracsStream);
    dBOrder = num2str(degree);
    wavetype = ['db' dBOrder];
    phi = zeros(1, signalLength);
    [phi_T, ~, ~] = wavefun(wavetype, 6);
    phi(1:length(phi_T))=phi_T;
    
    n_numVectors = 0:31;
    resolution = 64;
    
    % Create t
    % Construct the vectors of t^m, m = 0:degree-1 . 
    tVals = ones(signalLength, degree);
    for order = 2:degree
        tVals(:, order) = ((0:signalLength-1)/resolution).^(order - 1);
    end
    
    % Create shifted phi
    allPhi = zeros(length(phi),length(n_numVectors));
    for n = n_numVectors
        allPhi(:,n+1) = [zeros(1, n*resolution) phi(1:end - n*resolution)];
    end
    
    % Acquire coefficients
    coefficients =  (tVals.' * allPhi)./resolution;
    
    % Take xt and make yn
    yn =  x_diracsStream * allPhi;
    % Moments
    moments = coefficients * yn.';

end

