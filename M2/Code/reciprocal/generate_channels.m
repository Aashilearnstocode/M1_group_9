function [H_TX, H_RX] = generate_channels(p)
    % H_TX is R x N
    % H_RX is K x R
    
    % Rayleigh fading component
    H_TX_nlos = (randn(p.R, p.N) + 1i*randn(p.R, p.N))/sqrt(2);
    H_RX_nlos = (randn(p.K, p.R) + 1i*randn(p.K, p.R))/sqrt(2);
    
    % LoS component (random phase for simplicity)
    H_TX_los = exp(1i*2*pi*rand(p.R, p.N));
    H_RX_los = exp(1i*2*pi*rand(p.K, p.R));
    
    if p.Kf_dB == -Inf
        H_TX_small = H_TX_nlos;
        H_RX_small = H_RX_nlos;
    else
        H_TX_small = sqrt(p.Kf/(p.Kf+1))*H_TX_los + sqrt(1/(p.Kf+1))*H_TX_nlos;
        H_RX_small = sqrt(p.Kf/(p.Kf+1))*H_RX_los + sqrt(1/(p.Kf+1))*H_RX_nlos;
    end
    
    H_TX = sqrt(p.PL_BS_RIS) * H_TX_small;
    H_RX = sqrt(p.PL_RIS_user) * H_RX_small;
end
