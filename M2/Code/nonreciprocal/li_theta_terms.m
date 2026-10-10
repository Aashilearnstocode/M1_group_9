function [X, Y, Z] = li_theta_terms(G, Hr, V, iota, tau)
% LI_THETA_TERMS  Matrices of the BD-RIS sub-problem, Li et al. eq. (23),
%   reflective mode only (all users reflective, Phi_t = 0).
%   G  : R x N  (BS -> RIS)            = Aashi's H_TX (normalised)
%   Hr : K x R  (row k is h_k^H)       = Aashi's H_RX (normalised)
%   V  : N x K  precoder, tau/iota : K x 1
%   Sub-problem:  max_Theta  2 Re tr(Theta X) - tr(Theta Y Theta^H Z)
    taut = sqrt(1 + iota) .* tau;       % tilde-tau_k
    Gv   = G * V;                       % R x K, column k is g_k = G w_k
    X    = Gv * diag(conj(taut)) * Hr;  % sum_k conj(taut_k) g_k h_k^H
    Y    = Gv * Gv';                    % sum_p g_p g_p^H
    Z    = Hr' * diag(abs(tau).^2) * Hr;% sum_k |tau_k|^2 h_k h_k^H
end
