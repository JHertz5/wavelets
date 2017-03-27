function [ h_annihilatingFilter, tk_locations, ak_weights ] = annihilatingFilterMethodTLS( tau_moments )

    K = 2;
    N = length(tau_moments) - 1;
    % Construct S
    c = tau_moments(K+1:N);
    r = fliplr(tau_moments(1:K+1));
    S = toeplitz(c', r);
    % SVD to get V
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

