function [Theta, hist] = run_cga_optimizer(H_TX, H_RX, Theta_init, V, p, Rg, nu)
    max_iter = 8000;
    if isfield(p, 'eps_tol')
        eps_tol = p.eps_tol;
    else
        eps_tol = 1e-8;
    end
    
    Theta = Theta_init;
    
    % Initial variables
    [tau, y] = fp_variables(H_TX, H_RX, Theta, V, p);
    obj_val = objective(H_TX, H_RX, Theta, V, tau, y, p, nu);
    
    % Initial gradient and direction
    grad_eucl = gradient_theta(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu);
    r = tangent_projection(Theta, grad_eucl, Rg);
    
    % Initial direction is steepest ascent
    Xi = r; 
    
    hist.sumrate = zeros(max_iter, 1);
    hist.obj = zeros(max_iter, 1);
    hist.alpha = zeros(max_iter, 1);
    hist.gnorm = zeros(max_iter, 1);
    hist.armijo_max_hit = false(max_iter, 1);
    
    for iter = 1:max_iter
        % 1. Armijo line search
        [Theta_new, alpha, max_hit] = armijo_search(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu, obj_val, r, Xi);
        hist.alpha(iter) = alpha;
        hist.gnorm(iter) = norm(r, 'fro');
        hist.armijo_max_hit(iter) = max_hit;
        
        % 2. Update auxiliary variables
        [tau_new, y_new] = fp_variables(H_TX, H_RX, Theta_new, V, p);
        obj_new = objective(H_TX, H_RX, Theta_new, V, tau_new, y_new, p, nu);
        
        % Record history
        true_sum_rate = calculate_sumrate(tau_new);
        hist.sumrate(iter) = true_sum_rate;
        hist.obj(iter) = obj_new;
        
        % 3. Check convergence
        if abs(obj_new - obj_val) < eps_tol
            hist.sumrate = hist.sumrate(1:iter);
            hist.obj = hist.obj(1:iter);
            hist.alpha = hist.alpha(1:iter);
            hist.gnorm = hist.gnorm(1:iter);
            hist.armijo_max_hit = hist.armijo_max_hit(1:iter);
            Theta = Theta_new;
            break;
        end
        
        % 4. Compute new gradient and tangent projection
        grad_eucl_new = gradient_theta(H_TX, H_RX, Theta_new, V, tau_new, y_new, p, Rg, nu);
        r_new = tangent_projection(Theta_new, grad_eucl_new, Rg);
        
        % In manifold optimization, moving the previous gradient to the new tangent space 
        % requires vector transport. Projecting the direction when it is transported.
        Xi_transported = tangent_projection(Theta_new, Xi, Rg);
        r_transported = tangent_projection(Theta_new, r, Rg);
        
        % 5. Polak-Ribiere beta
        num_pr = real(trace(r_new' * (r_new - r_transported)));
        den_pr = real(trace(r_transported' * r_transported)) + 1e-16;
        beta = max(0, num_pr / den_pr);
        
        % 6. Ascent check
        Xi_new = r_new + beta * Xi_transported;
        
        % If <r_new, Xi_new> <= 0, reset to steepest ascent
        if real(trace(r_new' * Xi_new)) <= 0
            Xi_new = r_new;
        end
        
        % Update loop variables
        Theta = Theta_new;
        tau = tau_new;
        y = y_new;
        obj_val = obj_new;
        r = r_new;
        Xi = Xi_new;
    end
end
