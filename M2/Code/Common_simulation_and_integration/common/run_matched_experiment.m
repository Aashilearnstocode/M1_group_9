function res = run_matched_experiment(p)
% RUN_MATCHED_EXPERIMENT  Matched reciprocal-vs-general comparison (Member 4).
%
%   For every draw d = 1..p.num_draws
%       rng(p.seed_base + d)  ->  ONE channel (H_TX,H_RX), ONE start Theta_init, ONE V0
%   and for every architecture (SC, GC, FC; block size Rg = p.Rg_list)
%       reciprocal  (Member 2): run_cga_optimizer(nu = p.nu) + enforce_symmetry_unitarity
%       general-fv  (Member 3): general_fixedV_optimizer  (same CGA, nu = 0, V frozen)
%       general-jt  (Member 3): li_joint_optimizer        (Li Alg. 1, optional)
%   All final rates are re-evaluated HERE with common_sinr/common_sumrate, so
%   the numbers are guaranteed to come from one metric implementation.
%
%   Naming:  *_fv = rate with the frozen initial V0
%            *_v2 = rate after recomputing the MMSE V at the final Theta
%            *_jt = rate of Li's jointly optimised (Theta, V)
%   Like-for-like reciprocity gap:  gap_fv = gen_fv - rec_fv ,  gap_v2 = gen_v2 - rec_v2
%   (positive = reciprocity costs rate).  gen_jt is NOT comparable to rec_*
%   (it also optimises V) -- see M2/Code/common/README.md.
    setup_paths();
    nd = p.num_draws;   na = numel(p.Rg_list);
    nanm = @() nan(nd, na);

    res.p = p;
    res.rate_init = nan(nd, 1);
    res.chan_fro  = nan(nd, 2);                 % [||H_TX||_F  ||H_RX||_F] fingerprint
    res.rec_fv = nanm();  res.rec_v2 = nanm();  res.rec_it = nanm();  res.rec_time = nanm();
    res.rec_sym = nanm(); res.rec_uni = nanm();
    res.gen_fv = nanm();  res.gen_v2 = nanm();  res.gen_it = nanm();  res.gen_time = nanm();
    res.gen_sym = nanm(); res.gen_uni = nanm();
    res.gen_jt = nanm();  res.jt_it = nanm();   res.jt_time = nanm();
    res.rec_hist = cell(nd, na);  res.gen_hist = cell(nd, na);  res.jt_hist = cell(nd, na);
    res.errors = {};

    for d = 1:nd
        % ---- ONE shared channel / start / precoder per draw ----
        rng(p.seed_base + d, 'twister');
        [H_TX, H_RX] = common_generate_channels(p);
        Theta_init   = common_init_theta(p);
        V0           = common_mmse_beamformer(H_TX, H_RX, Theta_init, p);
        res.chan_fro(d, :) = [norm(H_TX, 'fro'), norm(H_RX, 'fro')];
        res.rate_init(d)   = common_sumrate(common_sinr(H_TX, H_RX, Theta_init, V0, p));

        for a = 1:na
            Rg = p.Rg_list(a);
            try
                % ---------- reciprocal (Member 2) ----------
                t = tic;
                [Th, h] = run_cga_optimizer(H_TX, H_RX, Theta_init, V0, p, Rg, p.nu);
                Th = enforce_symmetry_unitarity(Th, Rg);
                res.rec_time(d,a) = toc(t);
                res.rec_fv(d,a)  = common_sumrate(common_sinr(H_TX, H_RX, Th, V0, p));
                Vr = common_mmse_beamformer(H_TX, H_RX, Th, p);
                res.rec_v2(d,a)  = common_sumrate(common_sinr(H_TX, H_RX, Th, Vr, p));
                res.rec_it(d,a)  = numel(h.sumrate);
                res.rec_sym(d,a) = norm(Th - Th.', 'fro');
                res.rec_uni(d,a) = norm(Th'*Th - eye(p.R), 'fro');
                res.rec_hist{d,a} = h.sumrate;

                % ---------- general, V frozen (Member 3, matched) ----------
                t = tic;
                [Tg, hg] = general_fixedV_optimizer(H_TX, H_RX, Theta_init, V0, p, Rg);
                res.gen_time(d,a) = toc(t);
                res.gen_fv(d,a)  = common_sumrate(common_sinr(H_TX, H_RX, Tg, V0, p));
                Vg = common_mmse_beamformer(H_TX, H_RX, Tg, p);
                res.gen_v2(d,a)  = common_sumrate(common_sinr(H_TX, H_RX, Tg, Vg, p));
                res.gen_it(d,a)  = numel(hg.sumrate);
                res.gen_sym(d,a) = norm(Tg - Tg.', 'fro');
                res.gen_uni(d,a) = norm(Tg'*Tg - eye(p.R), 'fro');
                res.gen_hist{d,a} = hg.sumrate;

                % ---------- general, joint V/Theta (Member 3, paper-faithful) ----------
                if p.run_joint
                    t = tic;
                    [Tj, Vj, hj] = li_joint_optimizer(H_TX, H_RX, Theta_init, p, Rg);
                    res.jt_time(d,a) = toc(t);
                    res.gen_jt(d,a)  = common_sumrate(common_sinr(H_TX, H_RX, Tj, Vj, p));
                    res.jt_it(d,a)   = numel(hj.sumrate);
                    res.jt_hist{d,a} = hj.sumrate;
                end

                fprintf('draw %2d/%d %-2s (Rg=%2d) | rec %.3f (V2 %.3f) | gen-fv %.3f (V2 %.3f) | gap %.3f | gen-jt %.3f\n', ...
                    d, nd, p.arch_names{a}, Rg, res.rec_fv(d,a), res.rec_v2(d,a), ...
                    res.gen_fv(d,a), res.gen_v2(d,a), res.gen_fv(d,a) - res.rec_fv(d,a), res.gen_jt(d,a));
            catch err
                msg = sprintf('draw %d, %s: %s', d, p.arch_names{a}, err.message);
                res.errors{end+1} = msg;                       %#ok<AGROW>
                warning('run_matched_experiment:fail', '%s', msg);
            end
        end
    end

    % ---- like-for-like reciprocity gaps (positive = reciprocity costs rate) ----
    res.gap_fv     = res.gen_fv - res.rec_fv;
    res.gap_v2     = res.gen_v2 - res.rec_v2;
    res.gap_fv_pct = 100 * res.gap_fv ./ res.gen_fv;
    res.gap_v2_pct = 100 * res.gap_v2 ./ res.gen_v2;
    res.summary    = summarize_results(res);
end

function S = summarize_results(res)
% mean / standard error / n over valid draws, per architecture
    fields = {'rec_fv','rec_v2','gen_fv','gen_v2','gen_jt','gap_fv','gap_v2','gap_fv_pct','gap_v2_pct'};
    for f = 1:numel(fields)
        X = res.(fields{f});
        n = sum(~isnan(X), 1);
        S.(fields{f}).mean = mean(X, 1, 'omitnan');
        S.(fields{f}).se   = std(X, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
        S.(fields{f}).n    = n;
    end
end
