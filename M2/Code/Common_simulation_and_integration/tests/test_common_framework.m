% test_common_framework.m  --  sanity + parity checks for Member 4's common framework.
% Run:  run('M2/Code/tests/test_common_framework.m')   (every line prints PASS/FAIL)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));  setup_paths();
nfail = 0;
lbl = {'FAIL', 'PASS'};
chk = @(name, ok) fprintf('%s  %s\n', lbl{double(logical(ok)) + 1}, name);

p  = default_config();                 % common config (R = 32)
q  = config(32);                       % Member 2's config

%% 1. Config parity with reciprocal/config.m (every shared field)
shared = {'N','K','R','Pmax','N0','C0','d0','rho','d_BS_RIS','d_RIS_user','PL_BS_RIS','PL_RIS_user','Kf_dB','Kf'};
ok = true;
for i = 1:numel(shared), ok = ok && isequal(p.(shared{i}), q.(shared{i})); end
nfail = nfail + ~ok;  chk('default_config matches reciprocal/config.m (all shared fields)', ok);

%% 2. Channel parity: same seed -> identical channels AND identical RNG state afterwards
rng(7, 'twister');  [A, B] = common_generate_channels(p);  C1 = common_init_theta(p);
rng(7, 'twister');  [A2, B2] = generate_channels(q);       C2 = initialize_theta(q.R, 'SC', q.R);
ok = isequal(A, A2) && isequal(B, B2);
nfail = nfail + ~ok;  chk('common_generate_channels == generate_channels (same seed)', ok);
ok = isequal(C1, C2);
nfail = nfail + ~ok;  chk('common_init_theta == initialize_theta(R,''SC'',R) (RNG stream aligned)', ok);

%% 3. MMSE / SINR / sum-rate parity and power constraint
V  = common_mmse_beamformer(A, B, C1, p);
V2 = mmse_beamformer(A, B, C1, q);
ok = norm(V - V2, 'fro') < 1e-12 * norm(V2, 'fro');
nfail = nfail + ~ok;  chk('common_mmse_beamformer == mmse_beamformer', ok);
ok = abs(norm(V, 'fro')^2 - p.Pmax) < 1e-12 * p.Pmax;
nfail = nfail + ~ok;  chk('power constraint ||V||_F^2 = Pmax', ok);
s1 = common_sinr(A, B, C1, V, p);  s2 = calculate_sinr(A, B, C1, V, q);
ok = norm(s1 - s2) < 1e-10 * norm(s2);
nfail = nfail + ~ok;  chk('common_sinr == calculate_sinr', ok);
ok = abs(common_sumrate(s1) - calculate_sumrate(s2)) < 1e-10;
nfail = nfail + ~ok;  chk('common_sumrate == calculate_sumrate', ok);

%% 4. Reproducibility and independence of draws
rng(3, 'twister');  [X1, ~] = common_generate_channels(p);
rng(3, 'twister');  [X2, ~] = common_generate_channels(p);
rng(4, 'twister');  [X3, ~] = common_generate_channels(p);
ok = isequal(X1, X2) && ~isequal(X1, X3);
nfail = nfail + ~ok;  chk('same seed reproduces channel; different seed differs', ok);

%% 5. Path-loss scaling (large R so the sample mean is tight): E|h|^2 = PL, Rayleigh and Rician
for fad = {'rayleigh', 'rician'}
    pp = default_config(struct('R', 4000, 'fading', fad{1}));
    rng(1, 'twister');  [H1, H2] = common_generate_channels(pp);
    r1 = mean(abs(H1(:)).^2) / pp.PL_BS_RIS;   r2 = mean(abs(H2(:)).^2) / pp.PL_RIS_user;
    ok = abs(r1 - 1) < 0.05 && abs(r2 - 1) < 0.05;
    nfail = nfail + ~ok;  chk(sprintf('mean |h|^2 equals path loss (%s): ratios %.3f, %.3f', fad{1}, r1, r2), ok);
end

%% 6. Config overrides re-derive dependent fields
pp = default_config(struct('Pmax_dBm', 10, 'R', 16, 'G', 4));
ok = abs(pp.Pmax - 1e-2) < 1e-15 && isequal(pp.Rg_list, [1 4 16]);
nfail = nfail + ~ok;  chk('overrides re-derive Pmax and Rg_list', ok);

%% 7. Metric sanity: perfectly isolated users -> sum-rate = sum log2(1+|h|^2 P / N0)
pp = default_config();   pp.N = 2;   pp.K = 2;   pp.R = 2;
Vd = sqrt(0.05) * eye(2);
sn = common_sinr(eye(2), eye(2), eye(2), Vd, pp);
ok = norm(sn - 0.05 / pp.N0) < 1e-6 * (0.05 / pp.N0);
nfail = nfail + ~ok;  chk('SINR of an interference-free 2x2 case equals P/N0', ok);

if nfail == 0, fprintf('\nALL CHECKS PASSED\n'); else, fprintf('\nSOME CHECKS FAILED (%d)\n', nfail); end
