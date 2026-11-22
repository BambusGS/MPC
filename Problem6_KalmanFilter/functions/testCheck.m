function [nPass, nTest] = testCheck(ok, message, nPass, nTest)
    % Print PASS/FAIL for one check of testKalmanFilter.m and count it
    nTest = nTest + 1;
    if ok
        nPass = nPass + 1;
        fprintf('PASS  %s\n', message);
    else
        fprintf('FAIL  %s\n', message);
    end
end
