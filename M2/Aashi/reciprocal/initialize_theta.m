function Theta = initialize_theta(R, arch, G)
    % arch: 'SC', 'GC', 'FC'
    % G: number of groups (only used if GC)
    
    if strcmp(arch, 'SC')
        % Single connected: diagonal matrix with unit modulus
        phases = exp(1i*2*pi*rand(R, 1));
        Theta = diag(phases);
    elseif strcmp(arch, 'FC')
        % Fully connected: symmetric unitary
        [Q, ~] = qr((randn(R, R) + 1i*randn(R, R))/sqrt(2));
        Theta = Q * Q.';
    elseif strcmp(arch, 'GC')
        % Group connected
        Rg = R / G;
        Theta = zeros(R, R);
        for g = 1:G
            [Q, ~] = qr((randn(Rg, Rg) + 1i*randn(Rg, Rg))/sqrt(2));
            Theta((g-1)*Rg+1:g*Rg, (g-1)*Rg+1:g*Rg) = Q * Q.';
        end
    else
        error('Unknown architecture');
    end
end
