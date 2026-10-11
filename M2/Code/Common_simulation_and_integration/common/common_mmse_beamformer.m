function V = common_mmse_beamformer(H_TX, H_RX, Theta, p)
% COMMON_MMSE_BEAMFORMER  Regularised-ZF/MMSE precoder, scaled to ||V||_F^2 = Pmax.
%   V (N x K) = (E^H E + N0 I)^-1 E^H  with  E = H_RX * Theta * H_TX  (K x N).
    E = H_RX * Theta * H_TX;
    V_unscaled = (E' * E + p.N0 * eye(p.N)) \ E';
    V = sqrt(p.Pmax / norm(V_unscaled, 'fro')^2) * V_unscaled;
end
