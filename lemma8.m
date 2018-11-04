function H_final = lemma8(PI_theta, pi_t, A0, A1, A2, A3, B0, B1, B2, B3, K, V)
% Implementation of Lemma 8 from Robust Adaptive Tube Model Predictive
% Control

% Define important values

n_alpha = length(V(:,1));
p = 3;

vertices = compute_vertices(PI_theta,pi_t);
theta_hat_transpose = [ones(length(vertices(:,1)),1) vertices];

% set YALMIP decision variables
H = sdpvar(p + 1,n_alpha);

% set options for solver
options = sdpsettings('solver','gurobi');

for i = 1:n_alpha
    
    % set constraints
    Constraints = [];
    for j = 1:length(vertices(:,1))
        Constraints = [Constraints, theta_hat_transpose(j,:) * H >= 0];
    end

    Constraints = [Constraints, H * V == [V(i,:) * (A0 + B0 * K); V(i,:) * (A1 + B1 * K);
     V(i,:) * (A2 + B2 * K); V(i,:) * (A3 + B3 * K)]];

    % set objective that we want to minimize
    % Objective = max(theta_hat_transpose(j,:) * H * ones(length(H(1,:)),1));
    Objective = max(theta_hat_transpose * H * ones(length(H(1,:)),1));


    % solve the optimization
    sol = optimize(Constraints, Objective, options);

    H_final{i} = value(H);
    % H_final(i) = value(max(theta_hat_transpose(j,:) * H * ones(length(H(1,:)),1)));

end 