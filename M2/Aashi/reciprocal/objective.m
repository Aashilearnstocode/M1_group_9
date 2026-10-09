function obj = objective(H_TX, H_RX, Theta, V, tau, y, p, nu)
    % Evaluates the FP-transformed sum rate minus the symmetry penalty
    
    E = H_RX * Theta * H_TX; % K x N
    
    fp_sum_rate = 0;
    for k = 1:p.K
        e_k = E(k, :);
        
        sum_power = 0;
        for i = 1:p.K
            sum_power = sum_power + abs(e_k * V(:, i))^2;
        end
        
        term1 = log2(1 + tau(k));
        term2 = tau(k) / log(2);
        term3 = (1 + tau(k)) / log(2) * (2 * real(conj(y(k)) * (e_k * V(:, k))) - abs(y(k))^2 * (sum_power + p.N0));
        
        fp_sum_rate = fp_sum_rate + term1 - term2 + term3;
    end
    
    penalty = nu * norm(Theta - Theta.', 'fro')^2;
    
    obj = fp_sum_rate - penalty;
end
