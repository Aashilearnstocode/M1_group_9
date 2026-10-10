function [Theta, V, hist] = li_joint_optimizer(H_TX, H_RX, Theta_init, p, Rg, opts)
% LI_JOINT_OPTIMIZER  General (non-reciprocal-capable) BD-RIS baseline:
%   joint transmit-precoder / scattering-matrix design of
%   H. Li, S. Shen, B. Clerckx, "Beyond Diagonal Reconfigurable Intelligent
%   Surfaces: From Transmitting and Reflecting Modes to Single-, Group-, and
%   Fully-Connected Architectures" (Algorithm 1), reflective mode, K_t = 0.
%
%   Inputs  H_TX (R x N), H_RX (K x R)  -- exactly the outputs of Aashi's
%           generate_channels(p), so both methods see the same channel.
%           Theta_init : R x R block-unitary start (use the same one as Aashi)
%           p          : config struct (needs p.Pmax, p.N0)
%           Rg         : block size (1 = SC, R/G = GC, R = FC)
%           opts       : optional struct: max_outer (300), tol (1e-6 bps/Hz),
%                        max_cg (30), V_init (default: MMSE for Theta_init)
%   Outputs Theta (R x R, block-unitary, generally NOT symmetric),
%           V (N x K, original units, ||V||_F^2 <= p.Pmax),
%           hist.sumrate (per outer iteration, bps/Hz), hist.cg_iters.
%
%   Internally the channels are rescaled so that noise power = 1 and all
%   matrices are O(1); the SINR/sum-rate are invariant to this rescaling.
    if nargin < 6, opts = struct(); end
    max_outer = dflt(opts, 'max_outer', 300);
    tol       = dflt(opts, 'tol', 1e-6);
    max_cg    = dflt(opts, 'max_cg', 30);

    % ---- normalisation: noise -> 1, channels -> unit-variance scale ----
    aG = sqrt(mean(abs(H_TX(:)).^2));
    aH = sqrt(mean(abs(H_RX(:)).^2));
    G  = H_TX / aG;   Hr = H_RX / aH;
    scale = aG * aH / sqrt(p.N0);        % v_n = v * scale   (see derivation in README)
    Pn    = p.Pmax * scale^2;            % normalised power budget (= effective SNR budget)

    Theta = Theta_init;
    E = Hr * Theta * G;
    if isfield(opts, 'V_init')
        Vn = opts.V_init * scale;
    else
        Vu = (E' * E + eye(size(G, 2))) \ E';             % MMSE precoder, sigma^2 = 1
        Vn = sqrt(Pn) * Vu / norm(Vu, 'fro');
    end

    hist.sumrate = zeros(max_outer, 1);
    hist.cg_iters = zeros(max_outer, 1);
    rate_old = -inf;
    for it = 1:max_outer
        E = Hr * Theta * G;
        [iota, tau, sinr] = li_aux_update(E, Vn, 1);               % (17),(18)
        Vn = li_precoder_update(E, iota, tau, Pn);                 % (19),(20)
        [iota, tau] = li_aux_update(E, Vn, 1);                     % refresh with new V
        [X, Y, Z] = li_theta_terms(G, Hr, Vn, iota, tau);          % (23)
        [Theta, nit] = li_theta_update(Theta, X, Y, Z, Rg, max_cg);% Algorithm 2
        E = Hr * Theta * G;
        [~, ~, sinr] = li_aux_update(E, Vn, 1);
        rate = sum(log2(1 + sinr));
        hist.sumrate(it) = rate;
        hist.cg_iters(it) = nit;
        if abs(rate - rate_old) < tol, break; end
        rate_old = rate;
    end
    hist.sumrate  = hist.sumrate(1:it);
    hist.cg_iters = hist.cg_iters(1:it);
    Theta = project_unitary_blocks(Theta, Rg);                     % remove round-off
    V = Vn / scale;
end

function v = dflt(s, name, d)
    if isfield(s, name), v = s.(name); else, v = d; end
end
