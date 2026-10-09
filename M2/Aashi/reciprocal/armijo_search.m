function [Theta_new, alpha, max_hit] = armijo_search(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu, obj_val, T, Xi)
    % T is the Riemannian gradient (the tangent projection of Euclidean gradient).
    % Xi is the search direction (conjugate gradient direction).
    % We want to maximize, so f(Theta_new) >= f(Theta) + c * alpha * <T, Xi>
    
    max_steps = 200;
    max_hit = false;
    contraction = 0.75;
    c = 2e-11;
    
    % The paper scales the initial step to 1.
    % To prevent wild overshoot, we scale alpha0 by 1 / ||Xi||_F.
    alpha = 1 / norm(Xi, 'fro');
    if isinf(alpha) || isnan(alpha)
        alpha = 1;
    end
    
    % Directional derivative: real(trace(T' * Xi))
    % The Armijo condition for ascent: f(Theta + alpha*Xi) - f(Theta) >= c * alpha * <grad, Xi>
    inner_prod = real(trace(T' * Xi));
    
    for i = 1:max_steps
        Theta_new = retraction(Theta, Xi, alpha, Rg);
        obj_new = objective(H_TX, H_RX, Theta_new, V, tau, y, p, nu);
        
        if obj_new >= obj_val + c * alpha * inner_prod
            return;
        end
        
        alpha = alpha * contraction;
    end
    
    % If it fails to find an ascent step, return the last evaluated Theta_new
    max_hit = true;
end
