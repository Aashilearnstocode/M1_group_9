function f = li_theta_objective(Theta, X, Y, Z)
% Surrogate objective of the Theta sub-problem (to be MAXIMISED).
    f = 2 * real(trace(Theta * X)) - real(trace(Theta * Y * Theta' * Z));
end
