function [iota, tau, sinr] = li_aux_update(E, V, sig2)
% LI_AUX_UPDATE  Auxiliary variables of the fractional-programming transform
%   (Li et al., eqs. (17)-(18)).  E is the K x N effective channel H_RX*Theta*H_TX,
%   V is N x K, sig2 is the noise power in the same (normalised) units.
%     iota_k = gamma_k
%     tau_k  = sqrt(1+iota_k) * e_k v_k / (sum_p |e_k v_p|^2 + sig2)
    S    = E * V;                       % S(k,p) = e_k v_p
    P    = abs(S).^2;
    tot  = sum(P, 2) + sig2;
    sinr = diag(P) ./ (tot - diag(P));
    iota = sinr;
    tau  = sqrt(1 + iota) .* diag(S) ./ tot;
end
