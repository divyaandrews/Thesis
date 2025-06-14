function g = nogernd(alpha, beta)
    % Generate a uniform random number between 0 and 1
    u = rand();
    
    % Check if u is less than alpha
    if u < alpha
        % If true, return 0
        g = 0;
    else
        % Otherwise, generate a random number from the geometric distribution
        % with probability of success 1/beta plus 1
        p = 1 / beta;
        geomValue = geornd(p);
        g = 1 + geomValue;
    end
end