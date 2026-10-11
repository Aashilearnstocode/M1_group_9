function Theta = common_init_theta(p)
% COMMON_INIT_THETA  Shared starting scattering matrix: diagonal, unit-modulus,
%   random phase (single-connected). It is feasible for SC, GC and FC (and is
%   symmetric + unitary), so every method starts from the SAME point.
%   Draws R uniforms, same as reciprocal/initialize_theta(R,'SC',R).
    Theta = diag(exp(1i*2*pi*rand(p.R, 1)));
end
