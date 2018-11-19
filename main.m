function main()
tic % time the execution of the main function 

%% define important values for the problem
% length of prediction horizon 
N = 10;

A0 = [0.5 0.2;-0.1 0.6];
A1 = [0.042 0;0.072 0.03];
A2 = [0.015 0.019;0.009 0.035];
A3 = [0 0;0 0];
B0 = [0; 0.5];
B1 = [0;0];
B2 = [0;0];
B3 = [0.0397;0.059];

% initial condition:
x_t = [3; 6];

% next state of t
% x_t = [2.9; 3.2];

% initial control input
% u_t_1 = -1;  

% define sets for system parameters and distrubance
PI_theta = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
PI_w = [1 0;0 1;-1 0;0 -1];
pi_t = [1; 1; 1; 1; 1; 1];
pi_w = [0.1; 0.1; 0.1; 0.1];

% define F and G matrices to satisfy constraints on state and input
% F = [0 -3.33; 0 0];
% G = [0; 1];
F = [1/100 0; -1/100 0; 0 -10/3; 0 1/100; 0 0; 0 0];
G = [0; 0; 0; 0; 1; -1];


%% Offline Section of the Proposed Algorithm 
% Offline: given an initial parameter set estimate THETA_0, choose V and
% compute K, Y, H_c, H_1_hat, ..., H_n_alpha_hat, H_Q, and H_R

% Choose V matrix

% creates a V matrix with dimensions 2^5 x 2
V = generate_polytope3(2, 5);
% Take first 8 rows of V
V = V(1:8,:);
% append row of F that corresponds to row of 0s in G
V = [V; 0 -3.33]; % V is 9 x 2

% Calculate H_c and K, and lambda using lemma7():
[H_c, K, lambda] = lemma7(A0, A1, A2, A3, B0, B1, B2, B3, PI_theta, pi_t, F, G, V);


% Calculate H_1_hat, ..., H_n_alpha_hat by using lemma8():
H_hat = lemma8(PI_theta, pi_t, A0, A1, A2, A3, B0, B1, B2, B3, K, V);
% Access elements of H_hat by using H_hat{i}



% Calculate H_Q and H_R:
Q = [1 0;0 1];
R = 1;
[H_Q, H_R] = lemma10(V, Q, R, K);

% Offline section done! 

%% Online Section of the Proposed Algorithm 
state_evolution = figure;
state_sum = figure;

% Plot initial state
figure(state_evolution);
plot(x_t(1),x_t(2),'o','MarkerSize',5)
hold on 

% Initialise theta_hat_0 inTHETA_0
previous_point_estimate = [0 0 0]';

% Compute initial value of vertices for parameter set THETA_0
vertices = compute_vertices(PI_theta,pi_t);

mu = compute_mu(A1, A2, A3, B1, B2, B3);
% mu = 0.1;
% mu = 0.4;

% initial condition for point estimate
current_point_estimate = [0 0 0];

% create sum for states
sum_of_states = 0;

x = sdpvar(2,1);
for i = 1:500
    i
    sum_of_states = sum_of_states + norm(x_t)^2;
    figure(state_sum);
    hold on
    plot(i,sum_of_states,'x','MarkerSize',5)
    % Do the parameter set update with this function to get pi_t_plus_one
    if i ~= 1
        pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,optimal_control_input,x_t,PI_theta,PI_w,pi_t,pi_w);
        % calculate vertices of the newly updated parameter set:
        vertices = compute_vertices(PI_theta,(pi_t_plus_one)')
        % update the value of pi_t
        pi_t = pi_t_plus_one';
    end
    
    
    % update lambda_t:
    % lambda_t = update_lambda_t(vertices, H_hat);
    % disp(lambda_t)

    % store previous value for optimal_control_input
    if i ~= 1
        u_k_1 = optimal_control_input;
    end

    
    % Calculate point estimate
    if i ~= 1
        current_point_estimate = point_estimate(A0, A1, A2, A3, B0, B1, B2, B3, x_t, x_t_1, u_k_1, previous_point_estimate, PI_theta, pi_t, mu)
        previous_point_estimate = current_point_estimate;
    end
    
    
    
    % compute the optimal solution:
    theta_hat_transpose = [ones(length(vertices(:,1)),1) vertices];
    [optimal_cost, optimal_control_input, alpha_k_1] = compute_optimal_solution(A0, A1, A2, A3, B0, B1, B2, B3, N, H_c, G, theta_hat_transpose, H_hat, V, PI_w, pi_w, vertices, K, R, Q, x_t, current_point_estimate);

    
    % Store current value of state into the old value of state (update
    % value of old state)
    x_t_1 = x_t;
    % Update the state using the computed optimal control input and true
    % parameter value
    theta_used_to_update_state = [0.8 0.2 -0.5];
    [A_theta, B_theta] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,theta_used_to_update_state);
    w_t = [(randi([0 2000])-1000)/10000;(randi([0 2000])-1000)/10000];
    x_t = A_theta * x_t + B_theta * optimal_control_input + w_t;
    
    % Plot the newly computed state
    figure(state_evolution);
    plot(V*x<= alpha_k_1)
    plot(x_t(1),x_t(2),'o','MarkerSize',5)

    
    
end
% Title and labels for graph
title('Model Predictive Controller')
xlabel('x1') 
ylabel('x2') 

time_elapsed = toc
end 


