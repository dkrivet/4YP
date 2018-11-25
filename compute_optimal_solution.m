function [optimal_cost, optimal_control_input, alpha_k_1, predicted_v] = compute_optimal_solution(A0, A1, A2, A3, B0, B1, B2, B3, N, H_c, G, theta_hat_transpose, H_hat, V, PI_w, pi_w, vertices, K, R, Q, x_k, theta_hat)

m = length(vertices(:,1));
% use compute_w_bar to get the value of w_bar for use later on:
w_bar = compute_w_bar(PI_w, pi_w, V); % w_bar is array of length 9 

% set YALMIP decision variables:
% check the dimensions of these variables
v_k = sdpvar(N,1);
% the way this is defined, the columns of alpha_k make up alpha(0|k),
% alpha(1|k) , ... , alpha (N-1|k)
alpha_k = sdpvar(length(V(:,1)),N+1); % dimensions of 9 x 11




% need to figure out which value of x_k to use here, current value, or next
% value
[H, f_transpose] = construct_cost_matrices(A0, A1, A2, A3, B0, B1, B2, B3, K, N, R, Q, x_k, theta_hat);
% define objective here ...
Objective = v_k' * H * v_k + 2 * f_transpose * v_k;



% define constraints:
Constraints = [];
for i = 1:N
    Constraints = [Constraints, H_c * alpha_k(:,i) + G * v_k(i) <= ones(size(G))];
end

% ask about this line. Probably not right to just say alpha_k_plus_one is a
% vector of ones
% alpha_k_plus_one = ones(1,length(V(:,1)));
for i = 1:length(V(:,1))
    for j = 1:length(vertices(:,1))
        [A_theta_j, B_theta_j] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,vertices(j,:));
        for k = 1:N
            Constraints = [Constraints, theta_hat_transpose(j,:) * H_hat{i} * alpha_k(:,k) + V(i,:) * B_theta_j * v_k(k) + w_bar(i) <= alpha_k(i,k+1)];
        end 
    end
end

Constraints = [Constraints, alpha_k(:,1) >= V * x_k];
Constraints = [Constraints, H_c * alpha_k(:,N+1) <= ones(size(H_c * alpha_k(:,N+1)))];

for i = 1:length(V(:,1))
    for j = 1:m
        Constraints = [Constraints, alpha_k(i,N+1) >= theta_hat_transpose(j,:) * H_hat{i} * alpha_k(:,N+1) + w_bar(i)];
    end
end

% set solver options
options = sdpsettings('solver','gurobi','verbose',0);
% options = sdpsettings('solver','gurobi');

% solve the optimization
sol = optimize(Constraints, Objective, options);

optimal_cost = value(v_k' * H * v_k + 2 * f_transpose * v_k);
optimal_control_input = value(K * x_k + v_k(1));
% v_k(1)

    
% x = value(alpha_k(:,1));
% x(isnan(x)) = 0;
% alpha_k_0 = x;
alpha_k_1 = value(alpha_k(:,2));

predicted_v = value(v_k);

% value(H_c * alpha_k(:,10) + G * value(v_k(10,1)))

end 




