% test_gradient_fd.m
addpath(fileparts(mfilename('fullpath')));
rng(1, 'twister');

p = config(32);
nu = 1.0;
Rg_list = [1, 4, 32];
delta = 1e-6;

fprintf('--- Finite Difference Gradient Check ---\n');

[H_TX, H_RX] = generate_channels(p);

for Rg = Rg_list
    if Rg == 32
        arch = 'FC';
        G = 1;
    elseif Rg == 1
        arch = 'SC';
        G = p.R;
    else
        arch = 'GC';
        G = p.R / Rg;
    end
    
    % Create a non-symmetric block-unitary matrix to test the penalty gradient term
    Theta = zeros(p.R, p.R) + 1i * zeros(p.R, p.R);
    for g = 1:G
        idx = (g-1)*Rg + 1 : g*Rg;
        [Q, ~] = qr(randn(Rg, Rg) + 1i*randn(Rg, Rg));
        Theta(idx, idx) = Q;
    end
    
    V = mmse_beamformer(H_TX, H_RX, Theta, p);
    [tau, y] = fp_variables(H_TX, H_RX, Theta, V, p);
    
    % Compute analytic Euclidean gradient
    grad_analytic = gradient_theta(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu);
    
    grad_fd = zeros(p.R, p.R) + 1i * zeros(p.R, p.R);
    
    for r = 1:p.R
        g_idx = floor((r-1)/Rg) + 1;
        col_start = (g_idx - 1) * Rg + 1;
        col_end = g_idx * Rg;
        
        for c = col_start:col_end
            % Real part perturbation
            Theta_pos = Theta;
            Theta_pos(r, c) = Theta_pos(r, c) + delta;
            [tau_pos, y_pos] = fp_variables(H_TX, H_RX, Theta_pos, V, p);
            obj_pos = objective(H_TX, H_RX, Theta_pos, V, tau_pos, y_pos, p, nu);
            
            Theta_neg = Theta;
            Theta_neg(r, c) = Theta_neg(r, c) - delta;
            [tau_neg, y_neg] = fp_variables(H_TX, H_RX, Theta_neg, V, p);
            obj_neg = objective(H_TX, H_RX, Theta_neg, V, tau_neg, y_neg, p, nu);
            
            grad_fd(r, c) = grad_fd(r, c) + (obj_pos - obj_neg) / (2 * delta);
            
            % Imaginary part perturbation
            Theta_pos = Theta;
            Theta_pos(r, c) = Theta_pos(r, c) + 1i * delta;
            [tau_pos, y_pos] = fp_variables(H_TX, H_RX, Theta_pos, V, p);
            obj_pos = objective(H_TX, H_RX, Theta_pos, V, tau_pos, y_pos, p, nu);
            
            Theta_neg = Theta;
            Theta_neg(r, c) = Theta_neg(r, c) - 1i * delta;
            [tau_neg, y_neg] = fp_variables(H_TX, H_RX, Theta_neg, V, p);
            obj_neg = objective(H_TX, H_RX, Theta_neg, V, tau_neg, y_neg, p, nu);
            
            % Since it's a Wirtinger derivative w.r.t Theta*, the imaginary part contributes with +1i
            grad_fd(r, c) = grad_fd(r, c) + 1i * (obj_pos - obj_neg) / (2 * delta);
        end
    end
    
    err = norm(grad_analytic - grad_fd, 'fro') / norm(grad_fd, 'fro');
    fprintf('Rg = %2d: Relative Error = %e\n', Rg, err);
end
