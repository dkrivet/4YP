function [K_final, lambda] = calculate_K_from_V(A0, A1, A2, A3, B0, B1, B2, B3, PI_theta, pi_t, V)


n_alpha = length(V(:,1)); % # of rows in V
p = 3; % Because we have A0, A1, A2, A3


vertices = compute_vertices(PI_theta,pi_t);
theta_hat_transpose = [ones(length(vertices(:,1)),1) vertices];


% set YALMIP decision variables
lambda = sdpvar(1);
K = sdpvar(1,2);

H_hat = sdpvar((p + 1)*ones(1,n_alpha),n_alpha * ones(1,n_alpha));

% set constraints
Constraints = [];
for i = 1:n_alpha
    for j = 1: length(vertices(:,1))
        Constraints = [Constraints, lambda >= theta_hat_transpose(j,:) * H_hat{i} * ones(9,1), theta_hat_transpose(j,:) * H_hat{i} >= 0];
    end
end 

for i = 1:n_alpha
    Constraints = [Constraints, H_hat{i}*V==[V(i,:) * (A0 + B0 * K); V(i,:) * (A1 + B1 * K); V(i,:) * (A2 + B2 * K); V(i,:) * (A3 + B3 * K)]];
end
   
% set objective
Objective = lambda;

% set options for solver
options = sdpsettings('solver','gurobi');

% solve the optimization 
sol = optimize(Constraints, Objective, options);

K_final = value(K);

end

