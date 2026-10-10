% run_convergence_experiment.m
% Sum-rate vs iteration plot

addpath(fileparts(mfilename('fullpath')));

rng(1);

p = config(32); % R=32
p.N = 2;
p.K = 2;
p.Pmax_dBm = 20;
p.Pmax = 10^(p.Pmax_dBm/10) * 1e-3;
p.Kf_dB = -Inf; % Rayleigh
[H_TX, H_RX] = generate_channels(p);

% We can change this to 'FC' and 32 for the second run
arch = 'FC';
Rg = 32; % FC implies group size R
nu = 1.0;

Theta_init = initialize_theta(p.R, arch, 1); % G=1 for FC
V = mmse_beamformer(H_TX, H_RX, Theta_init, p);

initial_rate = calculate_sumrate(calculate_sinr(H_TX, H_RX, Theta_init, V, p));

fprintf('Running CGA Optimizer for %s architecture...\n', arch);
[Theta_opt, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);

% Prepend initial rate
hist.sumrate = [initial_rate; hist.sumrate];

Theta_final = enforce_symmetry_unitarity(Theta_opt, Rg);

sym_err = norm(Theta_final - Theta_final.', 'fro');
uni_err = norm(Theta_final' * Theta_final - eye(p.R), 'fro');

fprintf('Final Symmetry Error: %e\n', sym_err);
fprintf('Final Unitary Error: %e\n', uni_err);

% Evaluate final rate
[tau_final, ~] = fp_variables(H_TX, H_RX, Theta_final, V, p);
final_rate = calculate_sumrate(tau_final);
fprintf('Final Sum Rate: %.4f bps/Hz\n', final_rate);

% Diagnostic: Recompute V using Theta_final and check new rate
V2 = mmse_beamformer(H_TX, H_RX, Theta_final, p);
sinr2 = calculate_sinr(H_TX, H_RX, Theta_final, V2, p);
rate2 = calculate_sumrate(sinr2);
fprintf('Diagnostic - Sum Rate with recomputed V: %.4f bps/Hz\n', rate2);

% Plot and save
figure;
plot(0:length(hist.sumrate)-1, hist.sumrate, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Sum Rate (bps/Hz)');
title(sprintf('Convergence of Reciprocal BD-RIS (%s, R=32, 20 dBm, Rayleigh)', arch));
grid on;

saveas(gcf, sprintf('convergence_%s_R32.png', arch));
save(sprintf('convergence_results_%s.mat', arch), 'hist', 'Theta_final', 'sym_err', 'uni_err', 'p', 'final_rate');
