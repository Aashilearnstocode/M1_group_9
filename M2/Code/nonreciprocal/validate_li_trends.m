% validate_li_trends.m
% Qualitative validation against Li, Shen & Clerckx (Fig. 9 / Fig. 10 behaviour):
%   (i)  sum-rate increases with transmit power P,
%   (ii) fully-connected > group-connected > single-connected under Rayleigh fading,
%   (iii) the gap FC/GC over SC is large in Rayleigh fading.
% Paper setting (Fig. 9): N = K = 4, M = 32 cells, G = 8 groups, d_BI = 50 m, d_IU = 2.5 m,
% zeta_0 = -30 dB, eps = 2.2, noise -80 dBm.  The paper plots HYBRID-mode curves with K_t = K_r;
% here we run REFLECTIVE mode only (K_t = 0), so compare shape/ordering, not absolute numbers.
if ~exist('nd', 'var'), nd = 3; end                 % channel draws per point
P_dBm = [0 5 10 15 20];
add_reciprocal_path();
p = config(32);  p.N = 4;  p.K = 4;  p.Kf_dB = -Inf;       % Rayleigh
names = {'SC', 'GC (G=8)', 'FC'};  Rgs = [1 4 32];
opts = struct('max_outer', 150);
rate = zeros(numel(P_dBm), 3);
for ip = 1:numel(P_dBm)
    q = p;  q.Pmax_dBm = P_dBm(ip);  q.Pmax = 10^(P_dBm(ip)/10) * 1e-3;
    for d = 1:nd
        rng(d, 'twister');
        [H_TX, H_RX] = generate_channels(q);
        Th0 = initialize_theta(q.R, 'SC', q.R);
        for a = 1:3
            [Th, V] = li_joint_optimizer(H_TX, H_RX, Th0, q, Rgs(a), opts);
            rate(ip, a) = rate(ip, a) + calculate_sumrate(calculate_sinr(H_TX, H_RX, Th, V, q)) / nd;
        end
    end
    fprintf('P = %2d dBm | SC %.2f | GC %.2f | FC %.2f bps/Hz\n', P_dBm(ip), rate(ip, :));
end
if ~exist('results', 'dir'), mkdir('results'); end
save(fullfile('results', 'validate_li_trends.mat'), 'P_dBm', 'rate', 'names', 'nd');
try
    fig = figure('Visible', 'off');  plot(P_dBm, rate, '-o', 'LineWidth', 1.8);
    grid on; xlabel('Transmit power P (dBm)'); ylabel('Sum-rate (bps/Hz)');
    legend(names, 'Location', 'northwest');
    title('General BD-RIS baseline, reflective, Rayleigh (N=K=4, M=32)');
    saveas(fig, fullfile('results', 'validate_li_trends.png'));
catch err
    fprintf('Plotting skipped (%s).\n', err.message);
end
ok = all(all(diff(rate) > 0)) && all(rate(:,3) >= rate(:,2) - 1e-6) && all(rate(:,2) >= rate(:,1) - 1e-6);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('Trend check (increasing in P, FC >= GC >= SC): %s\n', tag);
