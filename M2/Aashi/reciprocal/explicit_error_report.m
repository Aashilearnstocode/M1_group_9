% explicit_error_report.m
addpath(fileparts(mfilename('fullpath')));

p = config(32);
arch = 'FC';
Rg = 32;
G = 1;
nu = 1.0;
p.eps_tol = 1e-8;

rng(1, 'twister');
[H_TX, H_RX] = generate_channels(p);

Theta_init = initialize_theta(p.R, arch, G);
V = mmse_beamformer(H_TX, H_RX, Theta_init, p);

fprintf('--- Error Report (FC, R=32) ---\n\n');

% 1. Initial Errors
err_sym_init = norm(Theta_init - Theta_init.', 'fro');
err_uni_init = norm(Theta_init' * Theta_init - eye(p.R), 'fro');
fprintf('Initial Theta:\n');
fprintf('  ||Theta - Theta^T||_F = %e\n', err_sym_init);
fprintf('  ||Theta^H * Theta - I||_F = %e\n\n', err_uni_init);

% 2. Optimize
[Theta_opt, ~] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);

% 3. Errors before projection
err_sym_opt = norm(Theta_opt - Theta_opt.', 'fro');
err_uni_opt = norm(Theta_opt' * Theta_opt - eye(p.R), 'fro');
fprintf('Optimized Theta (Before final projection):\n');
fprintf('  ||Theta - Theta^T||_F = %e\n', err_sym_opt);
fprintf('  ||Theta^H * Theta - I||_F = %e\n\n', err_uni_opt);

% 4. Errors after projection
Theta_final = enforce_symmetry_unitarity(Theta_opt, Rg);
err_sym_final = norm(Theta_final - Theta_final.', 'fro');
err_uni_final = norm(Theta_final' * Theta_final - eye(p.R), 'fro');
fprintf('Final Theta (After exact projection):\n');
fprintf('  ||Theta - Theta^T||_F = %e\n', err_sym_final);
fprintf('  ||Theta^H * Theta - I||_F = %e\n\n', err_uni_final);
