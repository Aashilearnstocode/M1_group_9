function T = tangent_projection(Theta, G, Rg)
    % Projects the Euclidean gradient G onto the tangent space of the unitary manifold
    
    R = size(Theta, 1);
    num_groups = R / Rg;
    
    T = zeros(R, R);
    
    for g = 1:num_groups
        idx = (g-1)*Rg + 1 : g*Rg;
        Theta_g = Theta(idx, idx);
        G_g = G(idx, idx);
        
        T_g = G_g - Theta_g * (Theta_g' * G_g + G_g' * Theta_g) / 2;
        
        T(idx, idx) = T_g;
    end
end
