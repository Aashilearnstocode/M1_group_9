function grad = gradient_theta(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu)
    % Computes the Euclidean gradient of the objective with respect to Theta
    
    grad = zeros(p.R, p.R);
    E = H_RX * Theta * H_TX;
    
    G = p.R / Rg;
    
    for g = 1:G
        idx = (g-1)*Rg + 1 : g*Rg;
        Theta_g = Theta(idx, idx);
        
        H_TX_g = H_TX(idx, :); % Rg x N
        H_RX_g = H_RX(:, idx); % K x Rg
        
        grad_sum_rate_g = zeros(Rg, Rg);
        
        for k = 1:p.K
            e_k = E(k, :);
            
            % Compute the summation over i
            sum_power_term = zeros(Rg, Rg);
            for i = 1:p.K
                hk_g = H_RX_g(k, :).'; % Rg x 1
                W_g_vi = H_TX_g * V(:, i); % Rg x 1
                
                term_i = conj(e_k * V(:, i)) * hk_g * (W_g_vi.'); % Rg x Rg
                sum_power_term = sum_power_term + term_i;
            end
            
            hk_g = H_RX_g(k, :).'; % Rg x 1
            W_g_vk = H_TX_g * V(:, k); % Rg x 1
            
            term1 = 2 * conj(y(k)) * hk_g * (W_g_vk.'); % Rg x Rg
            term2 = 2 * abs(y(k))^2 * sum_power_term; % Rg x Rg
            
            % Outer conjugation
            grad_k = ((1 + tau(k)) / log(2)) * conj(term1 - term2);
            
            grad_sum_rate_g = grad_sum_rate_g + grad_k;
        end
        
        grad_penalty_g = -4 * nu * (Theta_g - Theta_g.');
        
        grad(idx, idx) = grad_sum_rate_g + grad_penalty_g;
    end
end
