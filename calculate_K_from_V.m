function K_final = calculate_K_from_V(V)

% define some parameters of the problem:
A0 = [0.5 0.2;-0.1 0.6];
A1 = [0.042 0;0.072 0.03];
A2 = [0.015 0.019;0.009 0.035];
A3 = [0 0;0 0];
B0 = [0; 0.5];
B1 = [0;0];
B2 = [0;0];
B3 = [0.0397;0.059];

n_alpha = length(V(:,1)); % # of rows in V
p = 3; % Because we have A0, A1, A2, A3

PI_theta = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
pi_t = [1; 1; 1; 1; 1; 1];

vertices = con2vert(PI_theta,pi_t);
theta_hat_transpose = [ones(length(vertices(:,1)),1) vertices];


% set YALMIP decision variables
lambda = sdpvar(1);
K = sdpvar(1,2);
% H_1_hat = sdpvar(p + 1, n_alpha);
% H_2_hat = sdpvar(p + 1, n_alpha);
% H_3_hat = sdpvar(p + 1, n_alpha);
% H_4_hat = sdpvar(p + 1, n_alpha);
% H_5_hat = sdpvar(p + 1, n_alpha);
% H_6_hat = sdpvar(p + 1, n_alpha);
% H_7_hat = sdpvar(p + 1, n_alpha);
% H_8_hat = sdpvar(p + 1, n_alpha);
% H_9_hat = sdpvar(p + 1, n_alpha);
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

