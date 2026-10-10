function Theta = init_general_theta(R, Rg, mode)
% INIT_GENERAL_THETA  Initial scattering matrix for the general BD-RIS baseline.
%   mode = 'sc'   (default) diagonal, unit-modulus, random phase. This is the
%                 initialisation used by Li et al. (Sec. IV-F) and it is a valid
%                 starting point for SC, GC and FC alike.
%   mode = 'rand' random block-unitary (NOT symmetric) -- useful as a sanity
%                 test that the optimiser does not need a symmetric start.
    if nargin < 3, mode = 'sc'; end
    if strcmp(mode, 'sc')
        Theta = diag(exp(1i * 2*pi * rand(R, 1)));
    elseif strcmp(mode, 'rand')
        Theta = zeros(R, R);
        for g = 1:(R / Rg)
            idx = (g-1)*Rg + 1 : g*Rg;
            [Q, ~] = qr((randn(Rg) + 1i*randn(Rg)) / sqrt(2));
            Theta(idx, idx) = Q;
        end
    else
        error('init_general_theta: unknown mode "%s"', mode);
    end
end
