function [optimal_cost, optimal_control_input, alpha_k_current, alpha_k_to_return, predicted_v, sol] = compute_optimal_solution_V_form(current_time_step, A0, A1, A2, A3, B0, B1, B2, B3, N, F, G, V, PI_w, pi_w, vertices, K, R, Q, x_k, number_of_vertices, U_j, PI_theta, pi_theta, theta_hat, M0, previous_control_input, previous_v,h)

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
M = sdpvar(3,3,N-1);

Lambda = sdpvar(length(V(:,1)),6,N+1,number_of_vertices,'full');




% need to figure out which value of x_k to use here, current value, or next
% value
[H, f_transpose] = construct_cost_matrices(A0, A1, A2, A3, B0, B1, B2, B3, K, N, R, Q, x_k, theta_hat);
% define objective here ...
% if current_time_step <= 30
%     lambda = 10^6;
% else
%     lambda = 0;
% end

if mod(current_time_step,2) == 0
    lambda = 10^6;
else
    lambda = 0;
end
% lambda=10^6;

if current_time_step > 1
    Objective = v_k' * H * v_k + 2 * f_transpose * v_k - lambda * beta;
else
    Objective = v_k' * H * v_k + 2 * f_transpose * v_k;
end
% lambda = 10^4;
% Objective = v_k' * H * v_k + 2 * f_transpose * v_k - lambda * beta;  

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
        % Constraints = [Constraints, Lambda(:,:,i,j) * (pi_theta+h) <= alpha_k(:,i+1) - V * compute_little_d_of_x_and_u(A0, B0, U_j(:,:,j) * alpha_k(:,i), K * U_j(:,:,j) * alpha_k(:,i) + v_k(i)) - w_bar'];
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
    % Constraints = [Constraints, Lambda(:,:,N+1,j) * (pi_theta+h) <= alpha_k(:,N+1) - V * compute_little_d_of_x_and_u(A0, B0, U_j(:,:,j) * alpha_k(:,N+1), K * U_j(:,:,j) * alpha_k(:,N+1)) - w_bar'];
    Constraints = [Constraints, Lambda(:,:,N+1,j) * pi_theta <= alpha_k(:,N+1) - V * compute_little_d_of_x_and_u(A0, B0, U_j(:,:,j) * alpha_k(:,N+1), K * U_j(:,:,j) * alpha_k(:,N+1)) - w_bar'];
    % Constraint type: Element-wise inequality
    Constraints = [Constraints, Lambda(:,:,N+1,j) >= 0];
end



% D_x_u0 = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_k, previous_control_input);
% L_du = compute_L_of_du(B1, B2, B3, du);
% 
% % check this constraint is satisfied 
% % Constraint type: Matrix inequality
% Constraints = [Constraints, M0 + D_x_u0' * D_x_u0 + D_x_u0' * L_du + L_du' * D_x_u0 >= beta];
% % Constraint type: Element-wise inequality
% Constraints = [Constraints, beta >= 0];



% if current_time_step > 1 && current_time_step <=30
if current_time_step > 1 && mod(current_time_step,2) == 0
% if current_time_step > 1  
    % New constraints PE into future:
    [A_theta_hat, B_theta_hat] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,theta_hat);
    for i = 1:(N-1)
        v_hat(i) = previous_v(i+1);
    end
    
    x_hat = zeros(2,N);
    for i = 1:(N-1)
        if i == 1
            x_hat(:,i) = x_k;
        end
        u_hat(i) = K * x_hat(:,1) + v_hat(i);
        x_hat(:,i+1) = A_theta_hat * x_hat(:,i) + B_theta_hat * u_hat(i);
    end
    % v_hat, x_hat, and u_hat are now fully defined

%     for j = 1:length(U_j(:,:,1))    
%         sum = 0;
%         for i = 1:N-1
%             D_hat = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_hat(:,i), u_hat(i));
%             D_tilda = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, U_j(:,:,j)*alpha_k(:,i)-x_hat(i), u_hat(i));
%             sum = sum + D_hat'*D_hat + D_hat'*D_tilda + D_tilda'*D_hat;
%         end
%         Constraints = [Constraints, M0 + sum >= beta];
%     end
    for i =1:N-1
        sum = 0;
        for j = 1:length(U_j(:,:,1))
            D_hat = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_hat(:,i), u_hat(i));
            % D_tilda = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, U_j(:,:,j)*alpha_k(:,i)-x_hat(i), u_hat(i));
            D_tilda = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, U_j(:,:,j)*alpha_k(:,i)-x_hat(i), K*(U_j(:,:,j)*alpha_k(:,i)-x_hat(i)) + v_k(i)-v_hat(i));
            Constraints = [Constraints, D_hat'*D_hat + D_hat'*D_tilda + D_tilda'*D_hat >= M(:,:,i)];
        end
        sum = sum + M(:,:,i);
        Constraints = [Constraints, M(:,:,i) >= 0];
    end
    Constraints = [Constraints, sum >= beta];
    Constraints = [Constraints, beta>=0, beta<=100];
    %Constraints = [Constraints, beta<=100];
end


% set solver options
options = sdpsettings('solver','mosek','verbose',0, 'cachesolvers', 1);
% options = sdpsettings('solver','gurobi');

% solve the optimization
sol = optimize(Constraints, Objective, options);


optimal_cost = value(v_k' * H * v_k + 2 * f_transpose * v_k);
optimal_control_input = value(K * x_k + v_k(1));
alpha_k_current = value(alpha_k(:,1));
alpha_k_to_return = value(alpha_k);
predicted_v = value(v_k);





end 




