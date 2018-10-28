function H_c = lemma7()
% Set initial values

F = [0 -3.33; 0 0];
G = [0; 1];
% K = ...;

% use the generate_polytope3 script that Prof Cannon sent me
[V, theta, imax] = generate_polytope3(2, 5);
% this script creates a V matrix with dimensions (2^second_argument, 2)

% We want V with 9 rows, but we have to append a row from F to it so take
% first 8 rows of V
V = V(1:8,:);



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


