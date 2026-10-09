function sumrate = calculate_sumrate(sinr)
    sumrate = sum(log2(1 + sinr));
end
