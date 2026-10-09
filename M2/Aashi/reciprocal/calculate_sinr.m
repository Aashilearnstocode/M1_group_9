function sinr = calculate_sinr(H_TX, H_RX, Theta, V, p)
    % H_TX: R x N
    % H_RX: K x R
    % Theta: R x R
    % V: N x K
    
    E = H_RX * Theta * H_TX; % K x N
    
    sinr = zeros(p.K, 1);
    for k = 1:p.K
        h_k_tilde = E(k, :); % 1 x N
        signal_power = abs(h_k_tilde * V(:, k))^2;
        
        interf_power = 0;
        for j = 1:p.K
            if j ~= k
                interf_power = interf_power + abs(h_k_tilde * V(:, j))^2;
            end
        end
        
        sinr(k) = signal_power / (interf_power + p.N0);
    end
end
