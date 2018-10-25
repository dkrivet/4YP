% Implementation of Lemma 8 from Robust Adaptive Tube Model Predictive
% Control

% Define important values

PI_theta = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
PI_w = [1 0;0 1;-1 0;0 -1];
D = [(A1 * x_t_1+ B1* u_t_1) (A2 * x_t_1+ B2* u_t_1) (A3 * x_t_1+ B3* u_t_1)];
P_t = -PI_w * D;
vertices_of_theta = con2vert(PI_theta,pi_t);
% j will be chosen dynamically I assume
j = 2;
theta_hat_j = [ones(length(vertices_of_theta(1,:)),1) vertices_of_theta(j,:)'];
% we want theta_hat_j'

% set YALMIP decision variables
H = sdpvar(p + 1,n_alpha);

% set constraints
Constraints = [theta_hat_j * H, 1 <= j <= m, H * V = [V(i,:) * Phi_0; ----- ; V(i,:) * Phi_p]];

% set objective that we want to minimize
Objective = max(theta_hat_j * H * ones(length(H(1,:)),1));

% set options for solveer
options = sdpsettings('solver','gurobi');

% solve the optimization
sol = optimize(Constraints, Objective, options);