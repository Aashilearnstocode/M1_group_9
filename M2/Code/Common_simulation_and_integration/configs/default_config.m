function p = default_config(overrides)
% DEFAULT_CONFIG  Single source of truth for every simulation setting (Member 4).
%
%   p = default_config()               -> project M2 test configuration
%   p = default_config(struct('R',16,'Pmax_dBm',10))   -> override any field(s)
%
%   The returned struct deliberately uses the SAME field names as
%   reciprocal/config.m (p.N, p.K, p.R, p.Pmax, p.N0, p.PL_BS_RIS, p.PL_RIS_user,
%   p.Kf_dB, p.Kf, p.eps_tol), so Member 2's and Member 3's functions accept it
%   unchanged.  Derived fields (linear power, noise, path loss, Rg list) are
%   recomputed AFTER the overrides, so you only ever set the "_dBm"/"_dB" inputs.
%
%   ---- inputs --------------------------------------------------------------
%   N, K, R        TX antennas, single-antenna users, RIS elements
%   G              number of groups for group-connected (GC); Rg = R/G
%   Pmax_dBm       total BS transmit power budget, ||V||_F^2 <= Pmax   [dBm]
%   N0_dBm         noise power per user                                  [dBm]
%   C0_dB, d0, rho path loss: PL = C0 * (d/d0)^(-rho)
%   d_BS_RIS, d_RIS_user  link distances                                 [m]
%   fading         'rayleigh' or 'rician'   (rician uses Kf_dB, default 2 dB)
%   seed_base      draw d uses rng(seed_base + d, 'twister')
%   num_draws      number of matched channel draws
%   nu             reciprocal symmetry-penalty weight
%   eps_tol        CGA stopping tolerance (abs. change of surrogate objective)
%   run_joint      also run Li's joint V/Theta design (slower)
    if nargin < 1, overrides = struct(); end

    % ---------------- system ----------------
    p.N = 2;   p.K = 2;   p.R = 32;   p.G = 8;
    % ---------------- power / noise ----------------
    p.Pmax_dBm = 20;
    p.N0_dBm   = -80;
    % ---------------- path loss ----------------
    p.C0_dB = -30;   p.d0 = 1.0;   p.rho = 2.2;
    p.d_BS_RIS = 50; p.d_RIS_user = 2.5;
    % ---------------- fading ----------------
    p.fading = 'rayleigh';
    p.Kf_dB  = -Inf;
    % ---------------- experiment ----------------
    p.seed_base = 0;
    p.num_draws = 5;
    p.nu        = 1.0;
    p.eps_tol   = 1e-8;
    p.run_joint = true;

    % ---------------- apply overrides ----------------
    names = fieldnames(overrides);
    for i = 1:numel(names)
        p.(names{i}) = overrides.(names{i});
    end

    % ---------------- derived quantities ----------------
    switch lower(p.fading)
        case 'rayleigh'
            p.Kf_dB = -Inf;
        case 'rician'
            if isinf(p.Kf_dB), p.Kf_dB = 2; end       % default Rician K-factor
        otherwise
            error('default_config: fading must be ''rayleigh'' or ''rician''.');
    end
    if mod(p.R, p.G) ~= 0
        error('default_config: R (%d) must be divisible by G (%d).', p.R, p.G);
    end
    p.Pmax = 10^(p.Pmax_dBm/10) * 1e-3;                 % W
    p.N0   = 10^(p.N0_dBm/10)   * 1e-3;                 % W
    p.C0   = 10^(p.C0_dB/10);
    p.PL_BS_RIS   = common_pathloss(p.C0_dB, p.d_BS_RIS,   p.d0, p.rho);
    p.PL_RIS_user = common_pathloss(p.C0_dB, p.d_RIS_user, p.d0, p.rho);
    p.Kf = 10^(p.Kf_dB/10);                             % 0 for Rayleigh

    % ---------------- architectures (block size Rg) ----------------
    p.arch_names = {'SC', 'GC', 'FC'};
    p.Rg_list    = [1, p.R / p.G, p.R];                 % SC, GC, FC
end
