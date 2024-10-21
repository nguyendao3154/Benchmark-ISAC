cvx_begin
    variables x y
    % Ensure y is positive to maintain convexity
    x >= 10;
    y >= 2; 
    y<= 5;
    % Define the objective or constraint using x^2 / y
    expression f
    f = quad_over_lin(x,y);
    minimize(f)  % or use f in a constraint
cvx_end