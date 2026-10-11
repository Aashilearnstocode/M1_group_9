function R = common_sumrate(sinr)
% COMMON_SUMRATE  Sum-rate in bps/Hz:  sum_k log2(1 + gamma_k).
    R = sum(log2(1 + sinr));
end
