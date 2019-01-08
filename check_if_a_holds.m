function a_holds = check_if_a_holds(V, alpha)

% Define optimization variable
x = sdpvar(2,1);

for r = 1:length(V(:,1))
    V_to_mod = V;
    alpha_to_mod = alpha;
    V_to_mod(r,:) = [];
    alpha_to_mod(r) = [];
    
    % Define objective to minimize
    Objective = -V(r,:) * x;
    
    % Set constraints
    Constraints = [V_to_mod * x <= alpha_to_mod];
    
    % Set solver options
    options = sdpsettings('solver','gurobi','verbose',0);

    % Solve the optimization
    sol = optimize(Constraints, Objective, options);
    
    alpha_bar(r) = value(V(r,:) * x);
end

alpha_bar
alpha'
a_holds = all(alpha_bar' >= alpha);