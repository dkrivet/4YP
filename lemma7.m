function H_c = lemma7(F, G, K, V)
% Set initial values

% F = [0 -3.33; 0 0];
% G = [0; 1];
% K = ...
% V = ...    V has 9 rows


% Set YALMIP decision variables
H = sdpvar(1,n_alpha);


% set the value of n_alpha
% n_alpha = ... 

for i = 1:length(n_alpha)
    % Set constraints
    A = F + G * K;
    Constraints = [H >= 0, H * V == A(i,:)];
    
    % Set objective we want to minimize
    Objective = H * ones(length(H(1,:)), 1); % vector of ones must be the right size
    
    % set options for solver
    options = sdpsettings('solver','gurobi');

    % Solve the optimization
    sol = optimize(Constraints, Objective, options);

    
    H_c(i) = value(H * ones(length(H(1,:)), 1));
end
end


