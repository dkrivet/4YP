function optimal_cost = compute_optimal_solution(A0, A1, A2, A3, B0, B1, B2, B3, N, H_c, G, theta_hat_transpose, H_hat, V, PI_w, pi_w, vertices, K, R, Q, x_k, theta_hat)

m = length(vertices(:,1));
% use compute_w_bar to get the value of w_bar for use later on:
w_bar = compute_w_bar(PI_w, pi_w, V); % w_bar is array of length 9 

% set YALMIP decision variables:
% check the dimensions of these variables
v_k = sdpvar(N,1);
% the way this is defined, the columns of alpha_k make up alpha(0|k),
% alpha(1|k) , ... , alpha (N-1|k)
alpha_k = sdpvar(length(V(:,1)),N); % dimensions of 9 x 10
alpha_N_plus_one = sdpvar(length(V(:,1)), 1);



% need to figure out which value of x_k to use here, current value, or next
% value
[H, f_transpose] = construct_cost_matrices(A0, A1, A2, A3, B0, B1, B2, B3, K, N, R, Q, x_k, theta_hat);
% define objective here ...
Objective = v_k' * H * v_k + 2 * f_transpose * v_k;



% define constraints:
Constraints = [];
for i = 1:N
    Constraints = [Constraints, H_c * alpha_k(:,i) + G * v_k(i) <= ones(2,1)];
end

alpha_k_plus_one = ones(1,length(V(:,1)));
for i = 1:length(V(:,1))
    for j = 1:length(vertices(:,1))
        [A_theta_j, B_theta_j] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,vertices(j,:));
        Constraints = [Constraints, theta_hat_transpose(j,:) * H_hat{i} * alpha_k(:,i) + V(i,:) * B_theta_j * v_k(i) + w_bar(i) <= alpha_k_plus_one(i)];
    end
end

Constraints = [Constraints, alpha_k(1) >= V * x_k];
Constraints = [Constraints, H_c * alpha_N_plus_one <= ones(size(H_c * alpha_N_plus_one))];

for i = 1:length(V(:,1))
    for j = 1:m
        Constraints = [Constraints, alpha_N_plus_one(i) >= theta_hat_transpose(j,:) * H_hat{i} * alpha_N_plus_one + w_bar(i)];
    end
end

% set solver options
options = sdpsettings('solver','gurobi');

% solve the optimization
sol = optimize(Constraints, Objective, options);

optimal_cost = v_k' * H * v_k + 2 * f_transpose * v_k;

end 




