function [H_TX, H_RX] = common_generate_channels(p)
% COMMON_GENERATE_CHANNELS  BS->RIS and RIS->user channels (Member 4).
%   H_TX : R x N   (BS -> RIS)       H_RX : K x R   (RIS -> users)
%   Rayleigh if p.Kf_dB = -Inf, otherwise Rician with K-factor p.Kf and a
%   random-phase LoS component. Path loss is applied as sqrt(PL).
%
%   IMPORTANT (matched experiments): the random numbers are drawn in EXACTLY
%   the same order as reciprocal/generate_channels.m (nlos TX, nlos RX, los TX,
%   los RX -- the LoS draws are made even for Rayleigh), so
%       rng(d,'twister'); common_generate_channels(p)
%   returns the same channels as generate_channels(p) and leaves the generator
%   in the same state for the following initialize_theta call.
    H_TX_nlos = (randn(p.R, p.N) + 1i*randn(p.R, p.N)) / sqrt(2);
    H_RX_nlos = (randn(p.K, p.R) + 1i*randn(p.K, p.R)) / sqrt(2);
    H_TX_los  = exp(1i*2*pi*rand(p.R, p.N));
    H_RX_los  = exp(1i*2*pi*rand(p.K, p.R));

    if isinf(p.Kf_dB) && p.Kf_dB < 0                  % Rayleigh
        H_TX_small = H_TX_nlos;
        H_RX_small = H_RX_nlos;
    else                                              % Rician
        H_TX_small = sqrt(p.Kf/(p.Kf+1))*H_TX_los + sqrt(1/(p.Kf+1))*H_TX_nlos;
        H_RX_small = sqrt(p.Kf/(p.Kf+1))*H_RX_los + sqrt(1/(p.Kf+1))*H_RX_nlos;
    end

    H_TX = sqrt(p.PL_BS_RIS)   * H_TX_small;
    H_RX = sqrt(p.PL_RIS_user) * H_RX_small;
end
