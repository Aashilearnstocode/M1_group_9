function p = config(R_val)
    if nargin < 1
        R_val = 16;
    end
    
    p.N = 2; % TX antennas
    p.K = 2; % RX users
    p.R = R_val; % RIS elements
    
    p.Pmax_dBm = 20;
    p.Pmax = 10^(p.Pmax_dBm/10) * 1e-3; % Watts
    
    p.C0_dB = -30;
    p.C0 = 10^(p.C0_dB/10);
    p.d0 = 1.0;
    
    p.rho = 2.2;
    p.d_BS_RIS = 50;
    p.d_RIS_user = 2.5;
    
    p.N0_dBm = -80;
    p.N0 = 10^(p.N0_dBm/10) * 1e-3; % Watts
    
    % Path loss
    p.PL_BS_RIS = p.C0 * (p.d_BS_RIS / p.d0)^(-p.rho);
    p.PL_RIS_user = p.C0 * (p.d_RIS_user / p.d0)^(-p.rho);
    
    % Rician K-factor
    p.Kf_dB = -Inf; % Rayleigh by default, use 2 for Rician
    p.Kf = 10^(p.Kf_dB/10);
end
