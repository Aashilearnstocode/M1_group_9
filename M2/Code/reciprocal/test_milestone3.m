% test_milestone3.m
% Validates the complete CGA optimization loop

addpath(fileparts(mfilename('fullpath')));

rng(1);

p = config(16);
[H_TX, H_RX] = generate_channels(p);

arch = 'FC';
Rg = 16;
nu = 1.0;

Theta_init = initialize_theta(p.R, arch, 1);
V = mmse_beamformer(H_TX, H_RX, Theta_init, p);

initial_rate = calculate_sumrate(calculate_sinr(H_TX, H_RX, Theta_init, V, p));
fprintf('Initial Sum Rate: %.4f bps/Hz\n', initial_rate);

fprintf('Running CGA Optimizer for FC architecture...\n');
tic;
[Theta_opt, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);
toc;

% Pre-projection stats
pre_sym_err = norm(Theta_opt - Theta_opt.', 'fro');
fprintf('Pre-projection Symmetry Error: %e\n', pre_sym_err);
fprintf('Pre-projection Sum Rate: %.4f bps/Hz\n', hist.sumrate(end));

% Apply final symmetry/unitarity enforcement
Theta_final = enforce_symmetry_unitarity(Theta_opt, Rg);

% Final checks
sym_err = norm(Theta_final - Theta_final.', 'fro');
uni_err = norm(Theta_final' * Theta_final - eye(p.R), 'fro');

fprintf('\nFinal Symmetry Error: %e\n', sym_err);
fprintf('Final Unitary Error: %e\n', uni_err);
fprintf('Total Iterations: %d\n', length(hist.sumrate));

% Evaluate final rate
[tau_final, ~] = fp_variables(H_TX, H_RX, Theta_final, V, p);
final_rate = calculate_sumrate(tau_final);
fprintf('Final Sum Rate: %.4f bps/Hz\n', final_rate);
