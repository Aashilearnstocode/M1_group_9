function plot_matched_results(res, outdir)
% PLOT_MATCHED_RESULTS  Two standard figures for the matched comparison (Member 4).
%   matched_sumrate.png    mean +- SE sum-rate per architecture (rec vs general)
%   matched_convergence.png sum-rate vs iteration, draw 1, rec vs general (V frozen)
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    p = res.p;   S = res.summary;   na = numel(p.Rg_list);
    try
        % ---------- (1) bar chart ----------
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 400]);
        ax = axes(fig);  hold(ax, 'on');  grid(ax, 'on');
        M = [S.rec_fv.mean; S.gen_fv.mean; S.rec_v2.mean; S.gen_v2.mean]';
        E = [S.rec_fv.se;   S.gen_fv.se;   S.rec_v2.se;   S.gen_v2.se]';
        b = bar(ax, M, 'grouped');
        for i = 1:numel(b)                                  % error bars on each series
            x = b(i).XEndPoints;
            errorbar(ax, x, M(:, i)', E(:, i)', 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
        end
        set(ax, 'XTick', 1:na, 'XTickLabel', arrayfun(@(k) sprintf('%s (R_g=%d)', p.arch_names{k}, p.Rg_list(k)), 1:na, 'UniformOutput', false));
        ylabel(ax, 'Sum-rate (bps/Hz)');
        legend(ax, {'Reciprocal (V fixed)', 'General (V fixed)', 'Reciprocal (V recomputed)', 'General (V recomputed)'}, 'Location', 'northwest');
        title(ax, sprintf('Matched comparison: N=%d, K=%d, R=%d, P_{max}=%g dBm, %s, %d draws', ...
            p.N, p.K, p.R, p.Pmax_dBm, p.fading, p.num_draws));
        print(fig, fullfile(outdir, 'matched_sumrate.png'), '-dpng', '-r200');
        close(fig);

        % ---------- (2) convergence, draw 1 ----------
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 400]);
        ax = axes(fig);  hold(ax, 'on');  grid(ax, 'on');
        col = lines(na);
        for k = 1:na
            hr = res.rec_hist{1, k};  hg = res.gen_hist{1, k};
            if isempty(hr) || isempty(hg), continue; end
            plot(ax, 1:numel(hr), hr, '-',  'Color', col(k,:), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('%s reciprocal', p.arch_names{k}));
            plot(ax, 1:numel(hg), hg, '--', 'Color', col(k,:), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('%s general', p.arch_names{k}));
        end
        xlabel(ax, 'CGA iteration');  ylabel(ax, 'Sum-rate (bps/Hz)');
        legend(ax, 'Location', 'southeast');
        title(ax, sprintf('Convergence, draw 1 (V fixed): R=%d, P_{max}=%g dBm', p.R, p.Pmax_dBm));
        print(fig, fullfile(outdir, 'matched_convergence.png'), '-dpng', '-r200');
        close(fig);
    catch err
        warning('plot_matched_results:fail', 'Plotting skipped: %s', err.message);
    end
end
