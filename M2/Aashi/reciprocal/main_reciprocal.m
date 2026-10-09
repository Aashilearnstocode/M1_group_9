% % main_reciprocal.m
% Sanity check driver for Milestone 1

rng(1);

R_val = 16; % Default to 16, can change to 32
p = config(R_val);

fprintf('System Parameters: N=%d, K=%d, R=%d, Pmax=%.2f dBm\n', p.N, p.K, p.R, p.Pmax_dBm);

% Generate channels
[H_TX, H_RX] = generate_channels(p);

% Architectures to test
archs = {'SC', 'GC', 'FC'};
G_val = 2; % 2 groups for GC

for i = 1:length(archs)
    arch = archs{i};
    
    % Initialize Theta
    if strcmp(arch, 'GC')
        Theta = initialize_theta(p.R, arch, G_val);
        fprintf('\nArchitecture: %s (G=%d)\n', arch, G_val);
    else
        Theta = initialize_theta(p.R, arch, 1);
        fprintf('\nArchitecture: %s\n', arch);
    end
    
    % Compute MMSE beamformer
    V = mmse_beamformer(H_TX, H_RX, Theta, p);
    
    % Calculate SINR and sum rate
    sinr = calculate_sinr(H_TX, H_RX, Theta, V, p);
    sumrate = calculate_sumrate(sinr);
    
    fprintf('Sum Rate: %.4f bps/Hz\n', sumrate);
    
    % Symmetry check
    sym_err = norm(Theta - Theta.', 'fro');
    % Unitary check
    uni_err = norm(Theta' * Theta - eye(p.R), 'fro');
    
    fprintf('Symmetry Error: %e\n', sym_err);
    fprintf('Unitary Error: %e\n', uni_err);
end

%% ---- Self-checks (Milestone 1) ----
fprintf('\n---- Self-checks (Milestone 1) ----\n');
[H_TX, H_RX] = generate_channels(p);
Theta = initialize_theta(p.R, 'FC', 1);
V = mmse_beamformer(H_TX, H_RX, Theta, p);

% 1. Power constraint
fprintf('Power: %.6f W (Pmax = %.6f)\n', norm(V,'fro')^2, p.Pmax);

% 2. Independent SINR check (matrix form vs your loop)
E = H_RX*Theta*H_TX;  S = abs(E*V).^2;        % S(k,j) = |e_k v_j|^2
sinr2 = diag(S) ./ (sum(S,2) - diag(S) + p.N0);
fprintf('SINR mismatch: %e\n', norm(sinr2 - calculate_sinr(H_TX,H_RX,Theta,V,p)));

% 3. Theta must matter
r1 = calculate_sumrate(calculate_sinr(H_TX,H_RX,Theta,V,p));
Theta2 = initialize_theta(p.R,'FC',1);
r2 = calculate_sumrate(calculate_sinr(H_TX,H_RX,Theta2,V,p));
fprintf('Rate with Theta1: %.3f, with Theta2 (same V): %.3f\n', r1, r2);

% 4. Rate should grow with power
for dBm = [0 10 20]
    q = p; q.Pmax = 10^(dBm/10)*1e-3;
    Vq = mmse_beamformer(H_TX,H_RX,Theta,q);
    fprintf('Pmax %2d dBm -> %.3f bps/Hz\n', dBm, ...
        calculate_sumrate(calculate_sinr(H_TX,H_RX,Theta,Vq,q)));
end
