function alt_main(H_form, to_plot)
 % time the execution of the main function 

% Number of time steps to simulate:
sim_time = 30;
x = zeros(2,sim_time+1);
u = zeros(1, sim_time);

%% define important values for the problem
% length of prediction horizon 
N = 10;

% System matrices
A0 = [0.5 0.2;-0.1 0.6];
A1 = [0.042 0;0.072 0.03];
A2 = [0.015 0.019;0.009 0.035];
A3 = [0 0;0 0];
B0 = [0; 0.5];
B1 = [0;0];
B2 = [0;0];
B3 = [0.0397;0.059];

% initial condition:
x(:,1) = [3; 6];

% define sets for system parameters and disturbance
PI_theta = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
PI_w = [1 0;0 1;-1 0;0 -1];
pi_t = [1; 1; 1; 1; 1; 1];
pi_w = [0.1; 0.1; 0.1; 0.1];

% True theta value that we will try to converge to
true_theta = [0.8 0.2 -0.5];

% define F and G matrices to satisfy constraints on state and input
F = [1/10 0; -1/10 0; 0 -10/3; 0 1/10; 0 0; 0 0];
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
V = [V; 0 -10/3];

if ~H_form
    % Calculate U^j and R_j:
    alpha_bar = calculate_alpha_bar(V);
    indeces = rows_of_V_to_delete(alpha_bar);
    V(indeces,:) = [];
    tube_verts = con2vert(V, ones(length(V(:,1)),1));
    R_j = compute_R_j(V, tube_verts, ones(length(V(:,1)),1));
    U_j = calculate_U_j(R_j, V);
end
    
% Calculate stabilising gain K from our V matrix
[K, ~] = calculate_K_from_V(A0, A1, A2, A3, B0, B1, B2, B3, PI_theta, pi_t, V);

if H_form
    % Calculate H_c and K, and lambda using lemma7():
    H_c = lemma7(F, G, V, K);

    % Calculate H_1_hat, ..., H_n_alpha_hat by using lemma8():
    H_hat = lemma8(PI_theta, pi_t, A0, A1, A2, A3, B0, B1, B2, B3, K, V);
    % Access elements of H_hat by using H_hat{i}
end


% Calculate H_Q and H_R:
Q = [1 0;0 1];
R = 1;
% [H_Q, H_R] = lemma10(V, Q, R, K);

% Offline section done! 

%% Online Section of the Proposed Algorithm 
if to_plot
    state_evolution = figure;
    state_sum = figure;
    parameter_set = figure;

    % Plot initial state
    figure(state_evolution);
    plot(x(1,1),x(2,1),'o','MarkerSize',5)
    hold on 
end

% Initialise theta_hat_0 inTHETA_0
previous_point_estimate = [0 0 0]';

% Compute initial value of vertices for parameter set THETA_0
vertices = compute_vertices(PI_theta,pi_t);

% mu = compute_mu(A1, A2, A3, B1, B2, B3);
mu = 0.1;
% mu = 10;

% initial condition for point estimate
current_point_estimate = [0 0 0];

% create sum for states
if to_plot
    sum_of_states = 0;
    x_plot = sdpvar(2,1);
end

previous_v = 0;
for i = 1:sim_time
    tic
    i
    if to_plot
        sum_of_states = sum_of_states + norm(x(:,i))^2;
        figure(state_sum);
        hold on
        plot(i,sum_of_states,'x','MarkerSize',5)
    end
    
    % Do the parameter set update with this function to get pi_t_plus_one
    if i ~= 1
        pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x(:,i-1),u(i-1),x(:,i),PI_theta,PI_w,pi_t,pi_w);
        % calculate vertices of the newly updated parameter set:
        % remove semicolon at end of next line to output vertices 
        vertices = compute_vertices(PI_theta,(pi_t_plus_one)')
        % update the value of pi_t
        pi_t = pi_t_plus_one';
    end

    % compute size of parameter set    
    radial_size = compute_parameter_set_size(pi_t)
    
    % update lambda_t:
    % lambda_t = update_lambda_t(vertices, H_hat);
    % disp(lambda_t)


    
    % Calculate point estimate
    if i ~= 1
        % remove semicolon at end of line to output current point estimate
        current_point_estimate = point_estimate(A0, A1, A2, A3, B0, B1, B2, B3, x(:,i), x(:,i-1), u(i-1), previous_point_estimate, PI_theta, pi_t, mu);
        previous_point_estimate = current_point_estimate;
    end
    
    
    M0 = calculate_M0(x, u, i, A1, A2, A3, B1, B2, B3);
    % compute the optimal solution:
    theta_hat_transpose = [ones(length(vertices(:,1)),1) vertices];
    if i == 1
        previous_control_input = -1;
    else
        previous_control_input = u(i-1);
    end
    if H_form
        [~, u(i), ~, ~, ~, ~] = compute_optimal_solution(A0, A1, A2, A3, B0, B1, B2, B3, N, H_c, G, theta_hat_transpose, H_hat, V, PI_w, pi_w, vertices, K, R, Q, x(:,i), current_point_estimate);
    else
        [~, u(i), ~, alpha_k_1, predicted_v, ~] = compute_optimal_solution_V_form(i, A0, A1, A2, A3, B0, B1, B2, B3, N, F, G, V, PI_w, pi_w, vertices, K, R, Q, x(:,i), length(V(:,1)), U_j, PI_theta, pi_t, current_point_estimate, M0, previous_control_input,previous_v);
        previous_v = predicted_v;
    end
    
       
    % Store current value of state into the old value of state (update value of old state)
    % x_t_1 = x_t;
    
    % Update the state using the computed optimal control input and true parameter value
    [A_theta, B_theta] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,true_theta);
    % w_t = [(round(rand) - 0.5)/5;(round(rand) - 0.5)/5];
    w_t = [(randi([0 2000])-1000)/10000;(randi([0 2000])-1000)/10000];
    [~, idx] = max(w_t);
    w_t(idx) = (round(rand) - 0.5)/5;
    x(:,i+1) = A_theta * x(:,i) + B_theta * u(i) + w_t;

    
    % Plot the newly computed state
    % plot for all i|0 at initial time step 
    if to_plot
        figure(state_evolution);
        plot(V*x_plot<= alpha_k_1)
        plot(x(1,i+1),x(2,i+1),'o','MarkerSize',5)
    end
    
    time_elapsed = toc
    
end


if to_plot
    % Title and labels for graph
    title('Model Predictive Controller')
    xlabel('x1') 
    ylabel('x2') 

    % Plot terminal parameter set
    theta = sdpvar(3,1);
    figure(parameter_set);
    subplot(2,1,1);
    plot(PI_theta * theta <= ones(6,1))
    axis([-1 1 -1 1 -1 1])
    subplot(2,1,2);
    plot(PI_theta * theta <= pi_t)
    axis([-1 1 -1 1 -1 1])
end



end 