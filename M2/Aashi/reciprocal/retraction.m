function Theta_new = retraction(Theta, Xi, alpha, Rg)
    % Retraction onto the unitary manifold
    
    R = size(Theta, 1);
    num_groups = R / Rg;
    
    Theta_new = zeros(R, R);
    
    for g = 1:num_groups
        idx = (g-1)*Rg + 1 : g*Rg;
        Theta_g = Theta(idx, idx);
        Xi_g = Xi(idx, idx);
        
        Y = Theta_g + alpha * Xi_g;
        [Q, R_mat] = qr(Y);
        
        % Normalize Q to avoid arbitrary phase/sign flips from qr()
        D = diag(R_mat);
        ph = D ./ abs(D); % Phase of each diagonal element
        Q = Q * diag(ph);
        
        Theta_new(idx, idx) = Q;
    end
end
