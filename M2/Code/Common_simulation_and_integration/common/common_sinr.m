function sinr = common_sinr(H_TX, H_RX, Theta, V, p)
% COMMON_SINR  Per-user SINR (K x 1), linear.
%   gamma_k = |e_k v_k|^2 / ( sum_{j~=k} |e_k v_j|^2 + N0 ),  e_k = k-th row of H_RX*Theta*H_TX.
    E = H_RX * Theta * H_TX;          % K x N
    S = abs(E * V).^2;                % S(k,j) = |e_k v_j|^2
    sig = diag(S);
    sinr = sig ./ (sum(S, 2) - sig + p.N0);
end
