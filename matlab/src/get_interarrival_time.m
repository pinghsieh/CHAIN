function t = get_interarrival_time(lambda, number_of_mini_slots_per_unit_time)
% lambda = arrivals per unit time
% Poisson
    %t = exprnd(1/lambda);
% geometric
    p = lambda/number_of_mini_slots_per_unit_time;
    t = geornd(p)/number_of_mini_slots_per_unit_time;
end