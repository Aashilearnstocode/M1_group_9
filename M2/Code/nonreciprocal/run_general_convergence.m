% run_general_convergence.m
% Convergence experiment for the general (non-reciprocal-capable) BD-RIS baseline.
% Project test configuration: N=2, K=2, R=32, Pmax=20 dBm, Rayleigh, one channel draw.
% The call order  rng -> generate_channels -> initialize_theta -> mmse_beamformer
% is IDENTICAL to Aashi's run_monte_carlo.m / plot_convergence.m, so draw d here is
% the same channel + same start as draw d of her reciprocal runs.
add_reciprocal_path();
draw = 1;                                   % change to compare another matched draw
rng(draw, 'twister');

p = config(32);  p.N = 2;  p.K = 2;  p.Pmax_dBm = 20;
p.Pmax = 10^(p.Pmax_dBm/10) * 1e-3;  p.Kf_dB = -Inf;          % Rayleigh
[H_TX, H_RX] = generate_channels(p);
Theta_init = initialize_theta(p.R, 'SC', p.R);                 % shared diagonal start
V_init     = mmse_beamformer(H_TX, H_RX, Theta_init, p);       % shared initial precoder
rate_init  = calculate_sumrate(calculate_sinr(H_TX, H_RX, Theta_init, V_init, p));
fprintf('Initial sum-rate (SC start, MMSE V): %.4f bps/Hz\n', rate_init);

names = {'SC', 'GC (Rg=4, G=8)', 'FC'};   Rgs = [1 4 32];
res = struct();
for a = 1:3
    Rg = Rgs(a);
    % --- Mode A: paper-faithful joint design (Li et al., Algorithm 1) ---
    [Th_j, V_j, hj] = li_joint_optimizer(H_TX, H_RX, Theta_init, p, Rg);
    rj = calculate_sumrate(calculate_sinr(H_TX, H_RX, Th_j, V_j, p));
    % --- Mode B: matched, V frozen (Aashi's CGA with nu = 0) ---
    [Th_f, hf, of] = general_fixedV_optimizer(H_TX, H_RX, Theta_init, V_init, p, Rg);
    fprintf('%-15s joint: %.4f (%d outer its) | fixed-V: %.4f (%d its), V recomputed: %.4f | sym err %.2f\n', ...
        names{a}, rj, numel(hj.sumrate), of.rate_fixedV, numel(hf.sumrate), of.rate_recomputedV, of.sym_err);
    res(a).name = names{a};  res(a).Rg = Rg;
    res(a).joint_hist = [rate_init; hj.sumrate];   res(a).joint_rate = rj;
    res(a).fixedV_hist = [rate_init; hf.sumrate];  res(a).fixedV_out = of;
    res(a).Theta_joint = Th_j;  res(a).Theta_fixedV = Th_f;
end

if ~exist('results', 'dir'), mkdir('results'); end
save(fullfile('results', sprintf('general_convergence_draw%d.mat', draw)), 'res', 'p', 'draw', 'rate_init');

try
    fig = figure('Visible', 'off');
    subplot(1, 2, 1); hold on;
    for a = 1:3, plot(0:numel(res(a).joint_hist)-1, res(a).joint_hist, 'LineWidth', 1.8, 'DisplayName', res(a).name); end
    grid on; xlabel('Outer iteration (BCD)'); ylabel('Sum-rate (bps/Hz)');
    title('General BD-RIS, joint V/\Theta (Li Alg. 1)'); legend('Location', 'southeast');
    subplot(1, 2, 2); hold on;
    for a = 1:3, plot(0:numel(res(a).fixedV_hist)-1, res(a).fixedV_hist, 'LineWidth', 1.8, 'DisplayName', res(a).name); end
    grid on; xlabel('CGA iteration'); ylabel('Sum-rate (bps/Hz)');
    title('General BD-RIS, V fixed (matched to reciprocal)'); legend('Location', 'southeast');
    saveas(fig, fullfile('results', sprintf('general_convergence_draw%d.png', draw)));
catch err
    fprintf('Plotting skipped (%s). Results are in results/*.mat\n', err.message);
end
