% run_monte_carlo.m
addpath(fileparts(mfilename('fullpath')));

num_draws = 50;

% We should get p.Pmax_dBm and p.Kf_dB initialized via config
p = config(32);
p.N = 2;
p.K = 2;
p.Pmax_dBm = 20; 
p.Pmax = 10^(p.Pmax_dBm/10) * 1e-3;
p.Kf_dB = -Inf; % Rayleigh

arch = 'FC';
Rg = p.R;
G = 1;

nu = 1.0;

fprintf('Running %d Monte Carlo draws for %s, R=%d, Rg=%d...\n', num_draws, arch, p.R, Rg);

rates_pre = zeros(num_draws, 1);
rates = zeros(num_draws, 1);
rates_v2 = zeros(num_draws, 1);
iterations = zeros(num_draws, 1);
sym_errs = zeros(num_draws, 1);
uni_errs = zeros(num_draws, 1);
max_step_hits = zeros(num_draws, 1);

parfor i = 1:num_draws
    rng(i, 'twister');
    try
        [H_TX, H_RX] = generate_channels(p);
        
        % Fair canonical start
        Theta_init_SC = initialize_theta(p.R, 'SC', p.R);
        V = mmse_beamformer(H_TX, H_RX, Theta_init_SC, p);
        
        % Now init the architecture-specific Theta, but use the same V!
        % Or actually, the previous script used Theta_init_SC for all optimizers!
        % Yes, SC is a subset, so we can just start all of them at SC!
        Theta_init = Theta_init_SC;
        
        [Theta_opt, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);
        
        [tau_pre, ~] = fp_variables(H_TX, H_RX, Theta_opt, V, p);
        r_pre = calculate_sumrate(tau_pre);
        
        Theta_final = enforce_symmetry_unitarity(Theta_opt, Rg);
        
        [tau_final, ~] = fp_variables(H_TX, H_RX, Theta_final, V, p);
        r1 = calculate_sumrate(tau_final);
        
        V2 = mmse_beamformer(H_TX, H_RX, Theta_final, p);
        sinr2 = calculate_sinr(H_TX, H_RX, Theta_final, V2, p);
        r2 = calculate_sumrate(sinr2);
        
        sym_e = norm(Theta_final - Theta_final.', 'fro');
        uni_e = norm(Theta_final' * Theta_final - eye(p.R), 'fro');
        
        rates_pre(i) = r_pre;
        rates(i) = r1;
        rates_v2(i) = r2;
        iterations(i) = length(hist.sumrate);
        sym_errs(i) = sym_e;
        uni_errs(i) = uni_e;
        max_step_hits(i) = sum(hist.armijo_max_hit);
        
        fprintf('Draw %d: %d iters | Pre: %.4f | Post: %.4f | (V2: %.4f) | Armijo limit hits: %d\n', i, iterations(i), r_pre, r1, r2, max_step_hits(i));
    catch ME
        fprintf('Draw %d failed: %s\n', i, ME.message);
        rates(i) = NaN;
        rates_v2(i) = NaN;
        iterations(i) = NaN;
    end
end

% Compute statistics
valid = ~isnan(rates);
n_valid = sum(valid);

fprintf('\n=== Monte Carlo Results (%d valid draws) ===\n', n_valid);
if n_valid > 0
    fprintf('Iterations: Mean = %.1f, Std = %.1f\n', mean(iterations(valid)), std(iterations(valid)));
    fprintf('Sum Rate (Fixed V): Mean = %.4f, Std = %.4f, SE = %.4f, Min = %.4f, Max = %.4f\n', ...
        mean(rates(valid)), std(rates(valid)), std(rates(valid))/sqrt(n_valid), min(rates(valid)), max(rates(valid)));
    fprintf('Sum Rate (Recomputed V): Mean = %.4f, Std = %.4f, SE = %.4f, Min = %.4f, Max = %.4f\n', ...
        mean(rates_v2(valid)), std(rates_v2(valid)), std(rates_v2(valid))/sqrt(n_valid), min(rates_v2(valid)), max(rates_v2(valid)));
        
    diff = rates_v2(valid) - rates(valid);
    diff_proj = rates_pre(valid) - rates(valid);
    fprintf('Mean Rate Gain from Recomputing V: %.4f bps/Hz\n', mean(diff));
    fprintf('Mean Projection Cost: %e bps/Hz\n', mean(diff_proj));
    fprintf('Total Armijo 200-step limit hits: %d\n', sum(max_step_hits(valid)));
    fprintf('Max Symmetry Error: %e\n', max(sym_errs(valid)));
    fprintf('Max Unitary Error: %e\n', max(uni_errs(valid)));
end

save('mc_FC_R32.mat', 'rates_pre', 'rates', 'rates_v2', 'iterations', 'sym_errs', 'uni_errs', 'max_step_hits', 'p', 'arch', 'Rg');
