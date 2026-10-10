% plot_convergence.m
addpath(fileparts(mfilename('fullpath')));

p = config(32);
nu = 1.0;
eps_tol = 1e-8;
p.eps_tol = eps_tol;

architectures = {'SC', 'GC2', 'GC4', 'FC'};
% GC2 means Rg=2, GC4 means Rg=4
Rg_list = [1, 2, 4, 32];
num_draws = 20;

max_iters = 2000;
avg_histories = zeros(length(architectures), max_iters);
avg_rates = zeros(length(architectures), 1);
counts = zeros(length(architectures), max_iters);

fprintf('Running fair comparison over %d draws...\n', num_draws);

for d = 1:num_draws
    rng(d, 'twister');
    [H_TX, H_RX] = generate_channels(p);
    
    % Fair start: use SC init (diagonal) for all architectures
    Theta_init = initialize_theta(p.R, 'SC', p.R);
    V = mmse_beamformer(H_TX, H_RX, Theta_init, p);
    
    for i = 1:length(architectures)
        Rg = Rg_list(i);
        [Theta_opt, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);
        
        Theta_final = enforce_symmetry_unitarity(Theta_opt, Rg);
        [tau_final, ~] = fp_variables(H_TX, H_RX, Theta_final, V, p);
        rate = calculate_sumrate(tau_final);
        
        len = length(hist.sumrate);
        idx = min(len, max_iters);
        avg_histories(i, 1:idx) = avg_histories(i, 1:idx) + hist.sumrate(1:idx)';
        % Carry forward the last value
        if idx < max_iters
            avg_histories(i, idx+1:end) = avg_histories(i, idx+1:end) + hist.sumrate(end);
        end
        counts(i, :) = counts(i, :) + 1;
        avg_rates(i) = avg_rates(i) + rate;
    end
end

avg_histories = avg_histories ./ counts;
avg_rates = avg_rates / num_draws;

figure('Name', 'Convergence of BD-RIS Architectures (Fair Average)');
hold on;
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};
lines = {'-', '--', ':', '-.'};

for i = 1:length(architectures)
    % Find where it effectively flatlines to truncate plot gracefully
    effective_len = find(counts(i,:) > 0, 1, 'last');
    if isempty(effective_len), effective_len = max_iters; end
    % Just plot up to a reasonable x-axis (e.g. 1000)
    plot_len = min(effective_len, 1000);
    plot(1:plot_len, avg_histories(i, 1:plot_len), ...
        'Color', colors{i}, 'LineStyle', lines{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('%s (Avg Rate: %.2f)', architectures{i}, avg_rates(i)));
end

hold off;
grid on;
xlabel('Number of Iterations');
ylabel('Average Sum Rate (bps/Hz)');
title('Average Convergence of CGA for BD-RIS Architectures (R=32)');
legend('Location', 'southeast');
set(gca, 'FontSize', 12);

savefig('convergence_overlay.fig');
saveas(gcf, 'convergence_overlay.png');
fprintf('Saved averaged figure to convergence_overlay.png\n');
