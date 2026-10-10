function Gr = li_theta_gradient(Theta, X, Y, Z, Rg)
% Euclidean gradient of li_theta_objective (real-inner-product convention
% <A,B> = Re tr(A^H B)), kept only on the diagonal blocks of size Rg
% (the off-block entries are not free variables for SC/GC).
%   grad f = 2 ( X^H - Z Theta Y )          (cf. Li et al., eq. (30), sign flipped
%                                            because we maximise)
    Full = 2 * (X' - Z * Theta * Y);
    R  = size(Theta, 1);
    Gr = zeros(R, R);
    for g = 1:(R / Rg)
        idx = (g-1)*Rg + 1 : g*Rg;
        Gr(idx, idx) = Full(idx, idx);
    end
end
