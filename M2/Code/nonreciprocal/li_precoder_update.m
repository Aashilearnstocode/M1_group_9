function V = li_precoder_update(E, iota, tau, Pmax)
% LI_PRECODER_UPDATE  Transmit-precoder block of Li et al. (eqs. (19)-(20)).
%   w_k = ( sum_p hbar_p hbar_p^H + lambda I )^-1 * sqrt(1+iota_k) * hbar_k,
%   hbar_k = tau_k * htilde_k,  htilde_k = E(k,:)'.
%   lambda >= 0 is found by bisection so that ||W||_F^2 <= Pmax.
    [K, N] = size(E);
    Hb = E' * diag(tau);                         % N x K, column k = tau_k * htilde_k
    A  = Hb * Hb';                               % sum_p hbar_p hbar_p^H
    B  = Hb * diag(sqrt(1 + iota));              % N x K
    A  = (A + A') / 2;
    [U, D] = eig(A);  d = max(real(diag(D)), 0);
    C  = U' * B;
    pw = @(lam) sum(sum(abs(C).^2, 2) ./ (d + lam).^2);   % ||W(lambda)||_F^2
    if pw(0) <= Pmax && all(d > 0)
        lam = 0;
    else
        lo = 0;  hi = norm(B, 'fro') / sqrt(Pmax);        % pw(hi) <= Pmax
        for it = 1:200
            mid = (lo + hi) / 2;
            if pw(mid) > Pmax, lo = mid; else, hi = mid; end
            if (hi - lo) <= 1e-14 * max(hi, 1), break; end
        end
        lam = hi;
    end
    V = U * diag(1 ./ (d + lam)) * C;
end
