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

% We want to add to V a row that is a row of F that corresponds to the row
% of G that is 0

% Ive done it manually here but may want to do it dynamically later:
V = [V; 0 -3.33];
n_alpha = length(V(:,1));
K = calculate_K_from_V(V);



% Set YALMIP decision variables
H = sdpvar(1,n_alpha);


% set the value of n_alpha
n_alpha = length(V(:,1));


A = F + G * K;
for i = 1:length(A(:,1))
    % Set constraints
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


