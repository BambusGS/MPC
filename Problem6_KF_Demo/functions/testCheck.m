function [nPass, nTest] = testCheck(ok, message, nPass, nTest)
% TESTCHECK  Print PASS/FAIL for one check.
    nTest = nTest + 1;
    if ok
        nPass = nPass + 1;
        fprintf('PASS  %s\n', message);
    else
        fprintf('FAIL  %s\n', message);
    end
end
