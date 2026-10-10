function Theta_out = project_unitary_blocks(Theta, Rg)
% PROJECT_UNITARY_BLOCKS  Exact projection onto block-diagonal unitary matrices.
%   Non-reciprocal counterpart of Aashi's enforce_symmetry_unitarity.m:
%   it does the polar/SVD step but NO symmetrisation, so Theta is NOT forced
%   to satisfy Theta = Theta.'.
%   Rg = 1 -> SC (diagonal), Rg = R/G -> GC, Rg = R -> FC.
    R = size(Theta, 1);
    Theta_out = zeros(R, R);
    for g = 1:(R / Rg)
        idx = (g-1)*Rg + 1 : g*Rg;
        [U, ~, W] = svd(Theta(idx, idx));
        Theta_out(idx, idx) = U * W';
    end
end
