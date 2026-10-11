function PL = common_pathloss(C0_dB, d, d0, rho)
% COMMON_PATHLOSS  Large-scale power gain  PL = C0 * (d/d0)^(-rho)  (linear).
%   C0_dB : reference path loss at d0 [dB] (e.g. -30)
%   d     : link distance [m],  d0 : reference distance [m],  rho : exponent
    PL = 10^(C0_dB/10) * (d / d0)^(-rho);
end
