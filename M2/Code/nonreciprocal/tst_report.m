function tst_report(name, ok)
% Tiny helper for test_general_baseline.m: prints [PASS]/[FAIL] <name>.
    if ok, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %s\n', tag, name);
end
