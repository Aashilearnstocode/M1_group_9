function [Theta, nit] = li_theta_update(Theta, X, Y, Z, Rg, max_cg, c_armijo)
% LI_THETA_UPDATE  BD-RIS block of Li et al. (Algorithm 2): Riemannian
%   conjugate-gradient ascent on the (block-)unitary manifold.
%     * Riemannian gradient = tangent projection of the Euclidean gradient
%     * Polak-Ribiere beta (denominator uses the un-transported old gradient, as in Li)
%     * Armijo backtracking line search
%     * polar retraction  (Theta + a*Xi)(I + a^2 Xi^H Xi)^(-1/2)  ==  U*W' of its SVD
%   Rg = 1 -> SC, Rg = R/G -> GC, Rg = R -> FC.  NO symmetry constraint.
%   Simplification vs. the paper: all G blocks are updated jointly (product
%   manifold) rather than one group at a time; same fixed points, simpler code.
    if nargin < 6, max_cg = 30; end
    if nargin < 7, c_armijo = 1e-4; end
    f  = li_theta_objective(Theta, X, Y, Z);
    r  = tproj(Theta, li_theta_gradient(Theta, X, Y, Z, Rg), Rg);
    Xi = r;
    alpha = 1 / max(norm(Xi, 'fro'), eps);
    nit = 0;
    for it = 1:max_cg
        gn = norm(r, 'fro');
        if gn < 1e-12, break; end
        slope = real(trace(r' * Xi));
        % ---- Armijo backtracking ----
        ok = false;
        for ls = 1:50
            Tn = project_unitary_blocks(Theta + alpha * Xi, Rg);   % retraction
            fn = li_theta_objective(Tn, X, Y, Z);
            if fn >= f + c_armijo * alpha * slope, ok = true; break; end
            alpha = alpha / 2;
        end
        if ~ok, break; end
        nit = it;
        % ---- new gradient and Polak-Ribiere direction ----
        rn  = tproj(Tn, li_theta_gradient(Tn, X, Y, Z, Rg), Rg);
        rt  = tproj(Tn, r,  Rg);                    % transported old gradient
        Xt  = tproj(Tn, Xi, Rg);                    % transported old direction
        beta = max(0, real(trace(rn' * (rn - rt))) / (real(trace(r' * r)) + 1e-300));
        Xn  = rn + beta * Xt;
        if real(trace(rn' * Xn)) <= 0, Xn = rn; end % restart if not an ascent direction
        small = (fn - f) < 1e-13 * (1 + abs(f));
        Theta = Tn;  f = fn;  r = rn;  Xi = Xn;
        alpha = 2 * alpha;                          % be optimistic next time
        if small, break; end
    end
end

function T = tproj(Theta, Gam, Rg)
% Projection onto the tangent space of the (block) unitary manifold:
%   P(G) = G - Theta * sym(Theta^H G),  sym(A) = (A + A^H)/2   (block by block)
    R = size(Theta, 1);
    T = zeros(R, R);
    for g = 1:(R / Rg)
        idx = (g-1)*Rg + 1 : g*Rg;
        Tg = Theta(idx, idx);  Gg = Gam(idx, idx);
        M  = Tg' * Gg;
        T(idx, idx) = Gg - Tg * (M + M') / 2;
    end
end
