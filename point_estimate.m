function current_point_estimate = point_estimate(A0, A1, A2, A3, B0, B1, B2, B3, x_k, x_k_1, u_k_1, previous_point_estimate, PI_theta, pi_theta, mu)

D_x_1_u_1 = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_k_1, u_k_1);

[A_theta_hat_k_1, B_theta_hat_k_1] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,previous_point_estimate);
x_hat_1_at_k_1 = A_theta_hat_k_1 * x_k_1 + B_theta_hat_k_1 * u_k_1;

theta_k_tilda = previous_point_estimate + mu * D_x_1_u_1' * (x_k - x_hat_1_at_k_1);

% set YALMIP decision variables
theta = sdpvar(3,1);

% set YALMIP options
options = sdpsettings('solver','gurobi','verbose',0);

% set objective
Objective = (theta(1,1) - theta_k_tilda(1,1))^2 + (theta(2,1) - theta_k_tilda(2,1))^2 + (theta(3,1) - theta_k_tilda(3,1))^2;

% set constraints
Constraints = [PI_theta * theta <= pi_theta];

% solve the optimization
sol = optimize(Constraints, Objective, options);

current_point_estimate = value(theta);