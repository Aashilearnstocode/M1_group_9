function Theta_opt = enforce_symmetry_unitarity(Theta, Rg)
    % Enforces exact symmetry and unitarity at the end of the CGA optimization
    
    R = size(Theta, 1);
    num_groups = R / Rg;
    
    Theta_opt = zeros(R, R);
    
    for g = 1:num_groups
        idx = (g-1)*Rg + 1 : g*Rg;
        Theta_g = Theta(idx, idx);
        
        % Symmetrize
        Theta_sym = 0.5 * (Theta_g + Theta_g.');
        
        % Perform SVD
        [U, ~, V] = svd(Theta_sym);
        
        % Ensure unitarity
        Theta_opt(idx, idx) = U * V';
    end
end
