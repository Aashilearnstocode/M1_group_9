function V = mmse_beamformer(H_TX, H_RX, Theta, p)
    % Equivalent channel K x N
    E = H_RX * Theta * H_TX;
    
    % V = (E^H * E + N0 * I)^-1 * E^H
    V_unscaled = (E' * E + p.N0 * eye(p.N)) \ E';
    
    % Scale to ||V||_F^2 = Pmax
    scale = sqrt(p.Pmax / norm(V_unscaled, 'fro')^2);
    V = scale * V_unscaled;
end
