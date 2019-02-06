function [optimal_cost, optimal_control_input, alpha_k_current, alpha_k_1, predicted_v, sol] = compute_optimal_solution_V_form(A0, A1, A2, A3, B0, B1, B2, B3, N, F, G, V, PI_w, pi_w, vertices, K, R, Q, x_k, number_of_vertices, U_j, PI_theta, pi_theta, theta_hat, M0, previous_control_input)

m = length(vertices(:,1));
% use compute_w_bar to get the value of w_bar for use later on:
w_bar = compute_w_bar(PI_w, pi_w, V); % w_bar is array of length 9 

% set YALMIP decision variables:
% check the dimensions of these variables
v_k = sdpvar(N,1);
% the way this is defined, the columns of alpha_k make up alpha(0|k),
% alpha(1|k) , ... , alpha (N-1|k)
alpha_k = sdpvar(length(V(:,1)),N+1); % dimensions of 9 x 11

beta = sdpvar(1);
du = sdpvar(1);

Lambda = sdpvar(length(V(:,1)),6,N+1,number_of_vertices,'full');




% need to figure out which value of x_k to use here, current value, or next
% value
[H, f_transpose] = construct_cost_matrices(A0, A1, A2, A3, B0, B1, B2, B3, K, N, R, Q, x_k, theta_hat);
% define objective here ...
lambda = 10^6;
Objective = v_k' * H * v_k + 2 * f_transpose * v_k - lambda * beta;
% Objective = v_k' * H * v_k + 2 * f_transpose * v_k;


% define constraints:
Constraints = [];

% Constraint type: Element-wise inequality
for j = 1:number_of_vertices
    for i = 1:N
        Constraints = [Constraints, (F + G * K) * U_j(:,:,j) * alpha_k(:,i) + G * v_k(i) <= ones(size(G))];
    end
end



for i = 1:N
    for j = 1:number_of_vertices
        % Constraint type: Equality constraint
        Constraints = [Constraints, Lambda(:,:,i,j) * PI_theta == V * compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, U_j(:,:,j) * alpha_k(:,i), K * U_j(:,:,j) * alpha_k(:,i) + v_k(i))];
        % Constraint type: Element-wise inequality
        Constraints = [Constraints, Lambda(:,:,i,j) >= 0];
        % Constraint type: Element-wise inequality
        Constraints = [Constraints, Lambda(:,:,i,j) * pi_theta <= alpha_k(:,i+1) - V * compute_little_d_of_x_and_u(A0, B0, U_j(:,:,j) * alpha_k(:,i), K * U_j(:,:,j) * alpha_k(:,i) + v_k(i)) - w_bar'];
    end
end



% Constraint type: Element-wise inequality
Constraints = [Constraints, alpha_k(:,1) >= V * x_k];
% Constraints = [Constraints, H_c * alpha_k(:,N+1) <= ones(size(H_c * alpha_k(:,N+1)))];
for j = 1:number_of_vertices
    % Constraint type: Element-wise inequality
    Constraints = [Constraints, (F + G * K) * U_j(:,:,j) * alpha_k(:,N+1) <= 1];
end


for j = 1:number_of_vertices
    % Constraint type: Equality constraint
    Constraints = [Constraints, Lambda(:,:,N+1,j) * PI_theta == V * compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, U_j(:,:,j) * alpha_k(:,N+1), K * U_j(:,:,j) * alpha_k(:,N+1))];
    % Constraint type: Element-wise inequality
    Constraints = [Constraints, Lambda(:,:,N+1,j) * pi_theta <= alpha_k(:,N+1) - V * compute_little_d_of_x_and_u(A0, B0, U_j(:,:,j) * alpha_k(:,N+1), K * U_j(:,:,j) * alpha_k(:,N+1)) - w_bar'];
    % Constraint type: Element-wise inequality
    Constraints = [Constraints, Lambda(:,:,N+1,j) >= 0];
end



D_x_u0 = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_k, previous_control_input);
L_du = compute_L_of_du(B1, B2, B3, du);

% check this constraint is satisfied 
% Constraint type: Matrix inequality
Constraints = [Constraints, M0 + D_x_u0' * D_x_u0 + D_x_u0' * L_du + L_du' * D_x_u0 >= beta];
% Constraint type: Element-wise inequality
Constraints = [Constraints, beta >= 0];


% set solver options
options = sdpsettings('solver','mosek','verbose',0, 'cachesolvers', 1);
% options = sdpsettings('solver','gurobi');

% solve the optimization
sol = optimize(Constraints, Objective, options);

optimal_cost = value(v_k' * H * v_k + 2 * f_transpose * v_k);
optimal_control_input = value(K * x_k + v_k(1));
alpha_k_current = value(alpha_k(:,1));
alpha_k_1 = value(alpha_k(:,2));
predicted_v = value(v_k);



end 




