function optimal_cost = compute_optimal_solution(N, H_c, G, theta_hat_transpose, H_hat, V, PI_w, pi_w)

% use compute_w_bar to get the value of w_bar for use later on:
w_bar = compute_w_bar(PI_w, pi_w, V); % w_bar is array of length 9 

% set YALMIP decision variables:
% check the dimensions of these variables
v_k = sdpvar(N);
% the way this is defined, the columns of alpha_k make up alpha(0|k),
% alpha(1|k) , ... , alpha (N-1|k)
alpha_k = sdpvar(length(V(:,1)),N); % dimensions of 9 x 10
alpha_N_plus_one = sdpvar(length(V(:,1), 1);


% define objective here ...
Objective = v_k' * H * v_k + 2 * f' * v_k;



% define constraints:
Constraints = [];
for i = 1:N
    Constraints = [Constraints, H_c * alpha_k(i) + G * v_k(i) <= ones(2,1)];
end 

for i = 1:length(V(:,1))
    for j = 1:m
        Constraints = [Constraints, theta_hat_transpose(j,:) * H_hat{i} * alpha_k(i) + V(i,:) * B(theta(J)) * v_k(i) + w_bar(i) <= aplha_k_plus_one(i)];
    end
end

Constraints = [Constraints, alpha_k(1) >= V * x];
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

optimal_cost = v_k' * H * v_k + 2 * f' v_k;

end 




