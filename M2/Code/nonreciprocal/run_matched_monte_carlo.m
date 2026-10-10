% run_matched_monte_carlo.m
% Matched comparison reciprocal (Aashi) vs general (this folder) on identical draws.
%   rec_*  : reciprocal CGA, nu = 1, final symmetrise + polar projection (Aashi's pipeline)
%   gen_fv : general, V frozen,  same CGA machinery with nu = 0   <- like-for-like with rec_*
%   gen_jt : general, joint V/Theta (Li Algorithm 1)              <- paper-faithful baseline
% "fv" = rate with the frozen initial V;  "v2" = rate after recomputing the MMSE V at the end.
% The like-for-like reciprocity gap is  gen_fv - rec_fv  (and gen_v2 - rec_v2).
if ~exist('nd', 'var'), nd = 10; end            % number of channel draws (use 50 for the report)
add_reciprocal_path();
p = config(32);  p.N = 2;  p.K = 2;  p.Pmax_dBm = 20;
p.Pmax = 10^(p.Pmax_dBm/10) * 1e-3;  p.Kf_dB = -Inf;
names = {'SC', 'GC', 'FC'};  Rgs = [1 4 32];  nu = 1.0;
rec_fv = nan(nd,3); rec_v2 = nan(nd,3); gen_fv = nan(nd,3); gen_v2 = nan(nd,3); gen_jt = nan(nd,3);
rec_it = nan(nd,3); gen_it = nan(nd,3); jt_it = nan(nd,3);
for d = 1:nd
    rng(d, 'twister');
    [H_TX, H_RX] = generate_channels(p);
    Theta_init = initialize_theta(p.R, 'SC', p.R);
    V = mmse_beamformer(H_TX, H_RX, Theta_init, p);
    for a = 1:3
        Rg = Rgs(a);
        [Tr, hr] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);
        Tr = enforce_symmetry_unitarity(Tr, Rg);
        rec_fv(d,a) = calculate_sumrate(fp_variables(H_TX, H_RX, Tr, V, p));
        rec_v2(d,a) = calculate_sumrate(calculate_sinr(H_TX, H_RX, Tr, mmse_beamformer(H_TX, H_RX, Tr, p), p));
        rec_it(d,a) = numel(hr.sumrate);
        [~, hf, of] = general_fixedV_optimizer(H_TX, H_RX, Theta_init, V, p, Rg);
        gen_fv(d,a) = of.rate_fixedV;  gen_v2(d,a) = of.rate_recomputedV;  gen_it(d,a) = numel(hf.sumrate);
        [Tj, Vj, hj] = li_joint_optimizer(H_TX, H_RX, Theta_init, p, Rg);
        gen_jt(d,a) = calculate_sumrate(calculate_sinr(H_TX, H_RX, Tj, Vj, p));  jt_it(d,a) = numel(hj.sumrate);
        fprintf('draw %2d %-2s | rec %.3f (V2 %.3f) | gen fixed-V %.3f (V2 %.3f) | gen joint %.3f\n', ...
            d, names{a}, rec_fv(d,a), rec_v2(d,a), gen_fv(d,a), gen_v2(d,a), gen_jt(d,a));
    end
end
se = @(x) std(x) / sqrt(size(x,1));
fprintf('\n=== means over %d draws (bps/Hz)  [mean +- SE] ===\n', nd);
fprintf('%-4s | rec fixed-V   | gen fixed-V   | gap (gen-rec)  | gen joint\n', 'arch');
for a = 1:3
    gap = gen_fv(:,a) - rec_fv(:,a);
    fprintf('%-4s | %6.3f+-%.3f | %6.3f+-%.3f | %6.3f+-%.3f | %6.3f+-%.3f\n', names{a}, ...
        mean(rec_fv(:,a)), se(rec_fv(:,a)), mean(gen_fv(:,a)), se(gen_fv(:,a)), mean(gap), se(gap), mean(gen_jt(:,a)), se(gen_jt(:,a)));
end
if ~exist('results', 'dir'), mkdir('results'); end
save(fullfile('results', 'matched_mc.mat'), 'rec_fv', 'rec_v2', 'gen_fv', 'gen_v2', 'gen_jt', 'rec_it', 'gen_it', 'jt_it', 'names', 'Rgs', 'nd', 'p');
