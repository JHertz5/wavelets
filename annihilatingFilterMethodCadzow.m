function [ h_annihilatingFilter, tk_locations, ak_weights ] = annihilatingFilterMethodCadzow( tau_moments, iterations )
    
    K = 2;
    N = length(tau_moments) - 1;
    % Construct S
    c = tau_moments(K+1:N);
    r = fliplr(tau_moments(1:K+1));
    S = toeplitz(c', r);
    % SVD to get V
    for iteration = 1:iterations
        [U, D, V] = svd(S);
        % Remove smaller eigenvalues
        Dr = size(D, 1) - K;
        Dc = size(D, 2) - K;
        Ddash = [D(1:K, 1:K), zeros(K, Dc); zeros(Dr, K), zeros(Dr, Dc)];
        S = U * Ddash * V'; 
        % Make "Toeplitz"
        Pr = size(S, 1);
        Pc = size(S, 2);

        if Pr > Pc
            maxSize = Pr;
        else
            maxSize = Pc;
        end

        S_temp = zeros(maxSize, maxSize);
        for dIndex = -(Pr-1):1:Pc-1
            working = diag(S, dIndex);
            longth = maxSize - abs(dIndex);

            temp = mean(working);
            new = diag(temp * ones(longth, 1), dIndex);
            S_temp = S_temp + new;
        end
        S = S_temp(1:Pr, 1:Pc);
    end
    [~, ~, V] = svd(S);
    % H is last column of V
    h_annihilatingFilter = V(:, end);
    % As before
    tk_locations = roots(h_annihilatingFilter);
    
    % Solve the vandermonde system is:
    eqn_locationsMatrix = ones(K,K);
    for rowIndex = 2:K
        eqn_locationsMatrix(rowIndex,:) = tk_locations.' .^(rowIndex-1);
    end
    eqn_tauVect2 = tau_moments(1:K)';
    
    ak_weights = eqn_locationsMatrix \ eqn_tauVect2; % solve equation to find weights

end

