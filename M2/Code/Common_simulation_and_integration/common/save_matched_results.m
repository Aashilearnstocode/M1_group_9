function save_matched_results(res, outdir)
% SAVE_MATCHED_RESULTS  Standard result format for every matched experiment (Member 4).
%   outdir/matched_results.mat   full struct `res` (incl. histories, config)
%   outdir/per_draw.csv          one row per (draw, architecture)
%   outdir/summary.csv           mean/SE per architecture and metric
%   outdir/config_used.txt       every field of the config, human readable
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    p = res.p;   nd = p.num_draws;   na = numel(p.Rg_list);

    save(fullfile(outdir, 'matched_results.mat'), 'res', '-v7');

    % ---- per-draw table (long format) ----
    [D, A] = ndgrid(1:nd, 1:na);
    d = D(:);  a = A(:);
    T = table(d, p.seed_base + d, p.arch_names(a)', p.Rg_list(a)', ...
        res.chan_fro(d,1), res.chan_fro(d,2), res.rate_init(d), ...
        res.rec_fv(:), res.rec_v2(:), res.gen_fv(:), res.gen_v2(:), res.gen_jt(:), ...
        res.gap_fv(:), res.gap_v2(:), res.gap_fv_pct(:), ...
        res.rec_it(:), res.gen_it(:), res.jt_it(:), ...
        res.rec_sym(:), res.gen_sym(:), res.rec_uni(:), res.gen_uni(:), ...
        res.rec_time(:), res.gen_time(:), res.jt_time(:), ...
        'VariableNames', {'draw','seed','arch','Rg','H_TX_fro','H_RX_fro','rate_init', ...
        'rec_fv','rec_v2','gen_fv','gen_v2','gen_jt','gap_fv','gap_v2','gap_fv_pct', ...
        'rec_iters','gen_iters','jt_iters','rec_sym_err','gen_sym_err','rec_uni_err','gen_uni_err', ...
        'rec_time_s','gen_time_s','jt_time_s'});
    writetable(T, fullfile(outdir, 'per_draw.csv'));

    % ---- summary table ----
    metrics = fieldnames(res.summary);
    rows = {};
    for m = 1:numel(metrics)
        for k = 1:na
            s = res.summary.(metrics{m});
            rows(end+1, :) = {metrics{m}, p.arch_names{k}, p.Rg_list(k), s.mean(k), s.se(k), s.n(k)}; %#ok<AGROW>
        end
    end
    S = cell2table(rows, 'VariableNames', {'metric','arch','Rg','mean','se','n'});
    writetable(S, fullfile(outdir, 'summary.csv'));

    % ---- human-readable config ----
    fid = fopen(fullfile(outdir, 'config_used.txt'), 'w');
    names = fieldnames(p);
    for i = 1:numel(names)
        v = p.(names{i});
        if iscell(v), v = strjoin(v, ','); end
        if ischar(v), s = v; else, s = mat2str(v); end
        fprintf(fid, '%-12s = %s\n', names{i}, s);
    end
    fclose(fid);
    fprintf('Results saved to %s\n', outdir);
end
