% test_milestone2.m
% Validates Milestone 2 functions

rng(1);

p = config(16);
[H_TX, H_RX] = generate_channels(p);

arch = 'FC';
Rg = 16;
Theta = initialize_theta(p.R, arch, 1);
V = mmse_beamformer(H_TX, H_RX, Theta, p);

% 1. fp_variables
[tau, y] = fp_variables(H_TX, H_RX, Theta, V, p);

% 2. objective
nu = 1.0;
obj_val = objective(H_TX, H_RX, Theta, V, tau, y, p, nu);

% Verify obj_val matches the true sum rate calculation
true_sum_rate = calculate_sumrate(tau) - nu * norm(Theta - Theta.', 'fro')^2;
fprintf('Objective vs true sum rate diff: %e\n', abs(obj_val - true_sum_rate));

% 3. gradient_theta
grad = gradient_theta(H_TX, H_RX, Theta, V, tau, y, p, Rg, nu);

% Finite-difference test (Central difference)
epsilon = 1e-6;
grad_fd = zeros(p.R, p.R);
for i = 1:p.R
    for j = 1:p.R
        % Real part
        dT_re = zeros(p.R, p.R);
        dT_re(i, j) = epsilon;
        obj_plus_re = objective(H_TX, H_RX, Theta + dT_re, V, tau, y, p, nu);
        obj_minus_re = objective(H_TX, H_RX, Theta - dT_re, V, tau, y, p, nu);
        grad_re = (obj_plus_re - obj_minus_re) / (2 * epsilon);
        
        % Imaginary part
        dT_im = zeros(p.R, p.R);
        dT_im(i, j) = 1i * epsilon;
        obj_plus_im = objective(H_TX, H_RX, Theta + dT_im, V, tau, y, p, nu);
        obj_minus_im = objective(H_TX, H_RX, Theta - dT_im, V, tau, y, p, nu);
        grad_im = (obj_plus_im - obj_minus_im) / (2 * epsilon);
        
        % The paper's Euclidean gradient is 2 * Wirtinger derivative
        % Wirtinger derivative w.r.t Theta* is 0.5 * (d/dRe + i d/dIm)
        % Therefore Euclidean gradient = d/dRe + i d/dIm
        grad_fd(i, j) = grad_re + 1i * grad_im;
    end
end

rel_err = norm(grad - grad_fd, 'fro') / norm(grad, 'fro');
fprintf('Finite-difference relative error: %e\n', rel_err);

% 4. Tangent projection
T = tangent_projection(Theta, grad, Rg);
% The condition \Re{Theta^H T} = 0 in the paper means the Hermitian part is zero
tangent_check = norm(Theta' * T + (Theta' * T)', 'fro') / 2;
fprintf('Tangent check (Hermitian part of Theta^H T approx 0): %e\n', tangent_check);

% Check ascent direction
alpha_ascent = 1e-3;
Theta_ascent = retraction(Theta, T, alpha_ascent, Rg);
obj_ascent = objective(H_TX, H_RX, Theta_ascent, V, tau, y, p, nu);
fprintf('Objective change along T (alpha=1e-3): %e\n', obj_ascent - obj_val);

% 5. Retraction check
alpha = 0.1;
Theta_new = retraction(Theta, T, alpha, Rg);
retraction_check = norm(Theta_new * Theta_new' - eye(p.R), 'fro');
unitary_check_2 = norm(Theta_new' * Theta_new - eye(p.R), 'fro');
sym_retraction_check = norm(Theta_new - Theta_new.', 'fro');
fprintf('Retraction unitary check (Theta*Theta^H - I): %e\n', retraction_check);
fprintf('Retraction unitary check (Theta^H*Theta - I): %e\n', unitary_check_2);
fprintf('Retraction symmetry check: %e\n', sym_retraction_check);

% 6. Small-step consistency
alpha_tiny = 1e-10;
Theta_tiny = retraction(Theta, T, alpha_tiny, Rg);
consistency_check = norm(Theta_tiny - Theta, 'fro');
fprintf('Small-step consistency (Theta_new approx Theta): %e\n', consistency_check);
