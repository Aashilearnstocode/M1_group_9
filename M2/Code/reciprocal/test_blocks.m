% test_blocks.m
% Tests the M2 checks for Rg = 1, 4, 8

rng(1);
p = config(16);
[H_TX, H_RX] = generate_channels(p);

Rg_list = [1, 4, 8];
nu = 1.0;

for Rg = Rg_list
    fprintf('\n=== Testing Block Size Rg = %d ===\n', Rg);
    G = p.R / Rg;
    Theta = initialize_theta(p.R, 'GC', G);
    V = mmse_beamformer(H_TX, H_RX, Theta, p);
    
    [tau, y] = fp_variables(H_TX, H_RX, Theta, V, p);
    obj_val = objective(H_TX, H_RX, Theta, V, tau, y, p, nu);
    grad = gradient_theta(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu);
    
    % Finite difference
    epsilon = 1e-6;
    grad_fd = zeros(p.R, p.R);
    for g = 1:G
        idx = (g-1)*Rg + 1 : g*Rg;
        for i = 1:Rg
            for j = 1:Rg
                row = idx(i);
                col = idx(j);
                
                dT_re = zeros(p.R, p.R); dT_re(row, col) = epsilon;
                obj_plus_re = objective(H_TX, H_RX, Theta + dT_re, V, tau, y, p, nu);
                obj_minus_re = objective(H_TX, H_RX, Theta - dT_re, V, tau, y, p, nu);
                grad_re = (obj_plus_re - obj_minus_re) / (2 * epsilon);
                
                dT_im = zeros(p.R, p.R); dT_im(row, col) = 1i * epsilon;
                obj_plus_im = objective(H_TX, H_RX, Theta + dT_im, V, tau, y, p, nu);
                obj_minus_im = objective(H_TX, H_RX, Theta - dT_im, V, tau, y, p, nu);
                grad_im = (obj_plus_im - obj_minus_im) / (2 * epsilon);
                
                grad_fd(row, col) = grad_re + 1i * grad_im;
            end
        end
    end
    
    rel_err = norm(grad - grad_fd, 'fro') / norm(grad, 'fro');
    fprintf('Finite-difference relative error: %e\n', rel_err);
    
    T = tangent_projection(Theta, grad, Rg);
    tangent_check = norm(Theta' * T + (Theta' * T)', 'fro') / 2;
    fprintf('Tangent check (Hermitian part approx 0): %e\n', tangent_check);
end
