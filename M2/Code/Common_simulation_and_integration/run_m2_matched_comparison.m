% run_m2_matched_comparison.m  --  M2 matched reciprocal-vs-general comparison (Member 4)
%
% Usage (from anywhere):   run('M2/Code/run_m2_matched_comparison.m')
% Quick check (default):   num_draws = 3   (a few minutes)
% Report run:              num_draws = 50  (set BEFORE running, e.g. `num_draws = 50;`)
% Rician run:              fading = 'rician';
here = fileparts(mfilename('fullpath'));
addpath(here);  setup_paths();

if ~exist('num_draws', 'var'), num_draws = 3;           end
if ~exist('fading',    'var'), fading    = 'rayleigh';  end
if ~exist('run_joint', 'var'), run_joint = true;        end

p = default_config(struct('num_draws', num_draws, 'fading', fading, 'run_joint', run_joint));
fprintf('Matched experiment: N=%d K=%d R=%d G=%d Pmax=%g dBm %s, %d draws, seeds %d..%d\n', ...
    p.N, p.K, p.R, p.G, p.Pmax_dBm, p.fading, p.num_draws, p.seed_base+1, p.seed_base+p.num_draws);

res = run_matched_experiment(p);

outdir = fullfile(here, '..', 'Results', sprintf('matched_%s_R%d_%ddraws', p.fading, p.R, p.num_draws));
save_matched_results(res, outdir);
plot_matched_results(res, outdir);

% ---- console summary ----
S = res.summary;
fprintf('\n=== Means over %d draws (bps/Hz), mean +- SE ===\n', p.num_draws);
fprintf('%-4s | rec V-fixed   | gen V-fixed   | gap (gen-rec) | gap %%  | gen joint\n', 'arch');
for k = 1:numel(p.Rg_list)
    fprintf('%-4s | %5.3f+-%.3f | %5.3f+-%.3f | %5.3f+-%.3f | %5.2f%% | %5.3f\n', p.arch_names{k}, ...
        S.rec_fv.mean(k), S.rec_fv.se(k), S.gen_fv.mean(k), S.gen_fv.se(k), ...
        S.gap_fv.mean(k), S.gap_fv.se(k), S.gap_fv_pct.mean(k), S.gen_jt.mean(k));
end
if ~isempty(res.errors), fprintf('\n%d (draw, arch) runs FAILED:\n', numel(res.errors)); disp(res.errors'); end
