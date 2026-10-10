function [tau, y] = fp_variables(H_TX, H_RX, Theta, V, p)
    % Compute the auxiliary variables tau and y
    E = H_RX * Theta * H_TX; % K x N
    
    % tau is SINR
    tau = calculate_sinr(H_TX, H_RX, Theta, V, p);
    
    y = zeros(p.K, 1);
    for k = 1:p.K
        e_k = E(k, :); % 1 x N
        
        sum_power = 0;
        for i = 1:p.K
            sum_power = sum_power + abs(e_k * V(:, i))^2;
        end
        
        y(k) = (e_k * V(:, k)) / (sum_power + p.N0);
    end
end
