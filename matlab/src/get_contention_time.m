function t = get_contention_time(tau)
% Detrministic
    %t = tau;
% Exponential
    t = exprnd(tau);
end