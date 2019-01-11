function alpha_bar = calculate_alpha_bar(V)

% Define optimization variable
x = sdpvar(2,1);

alpha_bar = zeros(1, length(V(:,1)));

alpha = ones(1, length(V(:,1)));

for r = 1:length(V(:,1))
    % Create copies of V and alpha so I can remove rth row or element
    V_to_mod = V;
    alpha_to_mod = alpha;
    V_to_mod(r,:) = [];
    alpha_to_mod(r) = [];
    
    % Define objective to minimize
    Objective = -V(r,:) * x;
    
    % Set constraints
    Constraints = [V_to_mod * x <= alpha_to_mod'];
    
    % Set solver options
    options = sdpsettings('solver','gurobi','verbose',0);

    % Solve the optimization
    sol = optimize(Constraints, Objective, options);
    
    alpha_bar(r) = value(V(r,:) * x);
end

