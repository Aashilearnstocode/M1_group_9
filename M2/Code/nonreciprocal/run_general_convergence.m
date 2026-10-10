% run_general_convergence.m
% Convergence experiment for the general (non-reciprocal-capable) BD-RIS baseline.
% Project test configuration: N=2, K=2, R=32, Pmax=20 dBm, Rayleigh, one channel draw.
% The call order  rng -> generate_channels -> initialize_theta -> mmse_beamformer
% is IDENTICAL to reciprocal's run_monte_carlo.m / plot_convergence.m, so draw d here is
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

names = {'SC', 'GC (R_g=4, G=8)', 'FC'};   Rgs = [1 4 32];
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
    res(a).joint_hist  = [rate_init; hj.sumrate];   res(a).joint_rate = rj;
    res(a).fixedV_hist = [rate_init; hf.sumrate];   res(a).fixedV_out = of;
    res(a).Theta_joint = Th_j;  res(a).Theta_fixedV = Th_f;
end

if ~exist('results', 'dir'), mkdir('results'); end
save(fullfile('results', sprintf('general_convergence_draw%d.mat', draw)), 'res', 'p', 'draw', 'rate_init');

% ------------------------------------------------------------------------
% Plotting (to re-plot without rerunning the optimizers, put this block in
% its own script that starts with: load('results/general_convergence_draw1.mat'))
% ------------------------------------------------------------------------
try
    c  = lines(3);  lw = 1.8;
    fig = figure('Visible','off','Units','pixels','Position',[100 100 1000 400],'Color','w');
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

    % ---- (a) joint V/Theta ----
    ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
    for a = 1:3
        h = res(a).joint_hist;
        plot(ax1, 0:numel(h)-1, h, 'Color',c(a,:), 'LineWidth',lw, 'DisplayName',res(a).name);
    end
    xlabel(ax1, 'Outer iteration (BCD)');
    ylabel(ax1, 'Sum-rate (bps/Hz)');
    title(ax1, '(a) Joint V/\Theta (Li Alg. 1)');
    legend(ax1, 'Location','southeast', 'Box','off');
    xlim(ax1, [0 max(arrayfun(@(r) numel(r.joint_hist), res)) - 1]);

    % ---- (b) V fixed ----
    ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
    for a = 1:3
        h = res(a).fixedV_hist;  n = numel(h) - 1;
        plot(ax2, 0:n, h, 'Color',c(a,:), 'LineWidth',lw, 'DisplayName',res(a).name);
        plot(ax2, n, h(end), 'o', 'Color',c(a,:), 'MarkerFaceColor',c(a,:), ...
             'MarkerSize',4, 'HandleVisibility','off');      % mark where each run stopped
    end
    xlabel(ax2, 'CGA iteration');
    title(ax2, '(b) V fixed (matched to reciprocal)');
    legend(ax2, 'Location','southeast', 'Box','off');
    xlim(ax2, [0 max(arrayfun(@(r) numel(r.fixedV_hist), res)) - 1]);

    set([ax1 ax2], 'FontSize',12, 'GridAlpha',0.25, 'LineWidth',0.8);

    exportgraphics(fig, fullfile('results', sprintf('general_convergence_draw%d.png', draw)), 'Resolution',300);
    exportgraphics(fig, fullfile('results', sprintf('general_convergence_draw%d.pdf', draw)), 'ContentType','vector');
    close(fig);
catch err
    fprintf('Plotting skipped (%s). Results are in results/*.mat\n', err.message);
end