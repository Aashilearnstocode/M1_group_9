% test_general_baseline.m  --  sanity checks for the general (non-reciprocal-capable) BD-RIS baseline.
% Run from this folder.  Every check prints PASS/FAIL; a summary is printed at the end.
add_reciprocal_path();
nfail = 0;
chk = @(name, ok) tst_report(name, ok);

%% 1. Surrogate objective == direct per-user formula (Li eq. (21a))
rng(11);
p = config(8);  [H_TX, H_RX] = generate_channels(p);
G = H_TX / sqrt(mean(abs(H_TX(:)).^2));  Hr = H_RX / sqrt(mean(abs(H_RX(:)).^2));
Th = init_general_theta(8, 4, 'rand');
V  = (randn(2,2) + 1i*randn(2,2));  iota = rand(2,1)*3;  tau = randn(2,1) + 1i*randn(2,1);
[X, Y, Z] = li_theta_terms(G, Hr, V, iota, tau);
f1 = li_theta_objective(Th, X, Y, Z);
f2 = 0;
for k = 1:2
    taut = sqrt(1 + iota(k)) * tau(k);
    f2 = f2 + 2*real(conj(taut) * (Hr(k,:) * Th * G * V(:,k)));
    for q = 1:2
        f2 = f2 - abs(tau(k))^2 * abs(Hr(k,:) * Th * G * V(:,q))^2;
    end
end
ok = abs(f1 - f2) < 1e-10 * (1 + abs(f2));   nfail = nfail + ~ok;
chk(sprintf('sub-problem objective matches eq.(21a)  (|diff| = %.1e)', abs(f1-f2)), ok);

%% 2. Finite-difference check of the Euclidean gradient (blocks Rg = 1, 4, 8)
for Rg = [1 4 8]
    Th = init_general_theta(8, Rg, 'rand');
    Gr = li_theta_gradient(Th, X, Y, Z, Rg);
    D  = zeros(8);
    for g = 1:(8/Rg)
        idx = (g-1)*Rg+1 : g*Rg;
        D(idx, idx) = randn(Rg) + 1i*randn(Rg);
    end
    e = 1e-6;
    fd  = (li_theta_objective(Th + e*D, X, Y, Z) - li_theta_objective(Th - e*D, X, Y, Z)) / (2*e);
    an  = real(trace(Gr' * D));
    ok = abs(fd - an) < 1e-6 * (1 + abs(an));   nfail = nfail + ~ok;
    chk(sprintf('gradient vs finite difference, Rg=%d  (rel err = %.1e)', Rg, abs(fd-an)/(1+abs(an))), ok);
end

%% 3. Theta update: monotone ascent and unitarity
for Rg = [1 4 8]
    Th0 = init_general_theta(8, Rg, 'sc');
    f0  = li_theta_objective(Th0, X, Y, Z);
    Th1 = li_theta_update(Th0, X, Y, Z, Rg, 50);
    f1  = li_theta_objective(Th1, X, Y, Z);
    ue  = norm(Th1' * Th1 - eye(8), 'fro');
    ok = (f1 >= f0 - 1e-12) && ue < 1e-10;   nfail = nfail + ~ok;
    chk(sprintf('theta update ascends & stays unitary, Rg=%d  (f: %.4g -> %.4g, unitary err %.1e)', Rg, f0, f1, ue), ok);
end

%% 4. Joint optimiser on the project configuration (N=2, K=2, R=32, 20 dBm, Rayleigh)
rng(1);
p = config(32);  p.Kf_dB = -Inf;
[H_TX, H_RX] = generate_channels(p);
Th0 = initialize_theta(p.R, 'SC', p.R);
V0  = mmse_beamformer(H_TX, H_RX, Th0, p);
r0  = calculate_sumrate(calculate_sinr(H_TX, H_RX, Th0, V0, p));
res = zeros(3,1);  names = {'SC', 'GC(Rg=4)', 'FC'};  Rgs = [1 4 32];
for a = 1:3
    [Th, V, h] = li_joint_optimizer(H_TX, H_RX, Th0, p, Rgs(a));
    r_orig = calculate_sumrate(calculate_sinr(H_TX, H_RX, Th, V, p));   % Aashi's metric functions
    mono   = all(diff(h.sumrate) > -1e-9);
    pw     = norm(V, 'fro')^2 <= p.Pmax * (1 + 1e-9);
    ue     = norm(Th' * Th - eye(p.R), 'fro');
    % block structure respected?
    mask = kron(eye(p.R / Rgs(a)), ones(Rgs(a)));
    bs   = norm(Th .* (1 - mask), 'fro') < 1e-12;
    cons = abs(r_orig - h.sumrate(end)) < 1e-6;
    ok = mono && pw && ue < 1e-10 && bs && cons && r_orig > r0;   nfail = nfail + ~ok;
    chk(sprintf('%-9s joint: rate %.3f -> %.3f in %d outer its | monotone=%d power=%d unitary=%.0e blocks=%d metric-match=%d', ...
        names{a}, r0, r_orig, numel(h.sumrate), mono, pw, ue, bs, cons), ok);
    res(a) = r_orig;
end

%% 5. Matched fixed-V mode (Aashi's CGA with nu = 0): unitary, NOT forced symmetric
[Thf, hf, out] = general_fixedV_optimizer(H_TX, H_RX, Th0, V0, p, 32);
ok = out.uni_err < 1e-10 && out.rate_fixedV > r0;   nfail = nfail + ~ok;
chk(sprintf('fixed-V general FC: rate %.3f -> %.3f (V recomputed: %.3f), sym err %.2f (>0 expected), unitary err %.1e', ...
    r0, out.rate_fixedV, out.rate_recomputedV, out.sym_err, out.uni_err), ok);

%% 6. Qualitative paper behaviour (Li Fig. 9/10): FC >= GC >= SC on average (a few draws)
nd = 4;  rr = zeros(nd, 3);  fast = struct('max_outer', 120);
for d = 1:nd
    rng(100 + d);
    [H_TX, H_RX] = generate_channels(p);
    Th0 = initialize_theta(p.R, 'SC', p.R);
    for a = 1:3
        [Th, V] = li_joint_optimizer(H_TX, H_RX, Th0, p, Rgs(a), fast);
        rr(d, a) = calculate_sumrate(calculate_sinr(H_TX, H_RX, Th, V, p));
    end
end
m = mean(rr);
ok = m(3) >= m(2) - 1e-6 && m(2) >= m(1) - 1e-6;   nfail = nfail + ~ok;
chk(sprintf('mean over %d draws: SC %.3f <= GC %.3f <= FC %.3f bps/Hz', nd, m(1), m(2), m(3)), ok);

fprintf('\n%d check(s) failed.\n', nfail);
