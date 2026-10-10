function [Theta_final, hist, out] = general_fixedV_optimizer(H_TX, H_RX, Theta_init, V, p, Rg)
% GENERAL_FIXEDV_OPTIMIZER  "Matched" baseline: Aashi's reciprocal CGA machinery
%   (run_cga_optimizer: FP surrogate, V frozen, Armijo, retraction) with the
%   symmetry penalty switched off (nu = 0) and a plain unitary projection at the
%   end instead of enforce_symmetry_unitarity.  Same start, same V, same
%   tolerances as the reciprocal run  ==>  the only difference is the
%   reciprocity constraint, which is exactly what the project measures.
%   Requires reciprocal's functions on the path (see add_reciprocal_path.m).
    nu = 0;
    [Theta_opt, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu);
    Theta_final = project_unitary_blocks(Theta_opt, Rg);
    [tau, ~] = fp_variables(H_TX, H_RX, Theta_final, V, p);
    out.rate_fixedV = calculate_sumrate(tau);                       % V frozen
    V2 = mmse_beamformer(H_TX, H_RX, Theta_final, p);
    out.rate_recomputedV = calculate_sumrate(calculate_sinr(H_TX, H_RX, Theta_final, V2, p));
    out.sym_err = norm(Theta_final - Theta_final.', 'fro');         % expected > 0
    out.uni_err = norm(Theta_final' * Theta_final - eye(size(Theta_final,1)), 'fro');
end
