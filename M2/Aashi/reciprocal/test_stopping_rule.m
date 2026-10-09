% test_stopping_rule.m
addpath(fileparts(mfilename('fullpath')));

% We will test Draw 1 and Draw 29
draws = [1, 29];
tols = [1e-8, 1e-10, 1e-12];

p = config(32);
arch = 'FC';
Rg = 32;
nu = 1.0;

for d = draws
    fprintf('\n=== Testing Draw %d ===\n', d);
    % Always use 'twister' for reproducibility
    rng(d, 'twister');
    [H_TX, H_RX] = generate_channels(p);
    Theta_init = initialize_theta(p.R, arch, 1);
    V = mmse_beamformer(H_TX, H_RX, Theta_init, p);
    
    for tol = tols
        p.eps_tol = tol;
        
        tic;
        [Theta_opt, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);
        t_elapsed = toc;
        
        Theta_final = enforce_symmetry_unitarity(Theta_opt, Rg);
        [tau_final, ~] = fp_variables(H_TX, H_RX, Theta_final, V, p);
        rate = calculate_sumrate(tau_final);
        
        fprintf('Tol %e: %d iters | Rate = %.4f | Time = %.2fs | alpha = %e | gnorm = %e\n', tol, length(hist.sumrate), rate, t_elapsed, hist.alpha(end), hist.gnorm(end));
    end
end
