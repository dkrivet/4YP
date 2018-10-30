function main()
%% define important values for the problem
A0 = [0.5 0.2;-0.1 0.6];
A1 = [0.042 0;0.072 0.03];
A2 = [0.015 0.019;0.009 0.035];
A3 = [0 0;0 0];
B0 = [0; 0.5];
B1 = [0;0];
B2 = [0;0];
B3 = [0.0397;0.059];

% initial condition:
x_t_1 = [3; 6];

% next state of t
x_t = [2.9; 3.2];

% initial control input
u_t_1 = -1;  

% define sets for system parameters and distrubance
PI_theta = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
PI_w = [1 0;0 1;-1 0;0 -1];
pi_t = [1; 1; 1; 1; 1; 1];
pi_w = [0.1; 0.1; 0.1; 0.1];

% define F and G matrices to satisfy constraints on state and input
F = [0 -3.33; 0 0];
G = [0; 1];

%% Offline Section of the Proposed Algorithm 
% Offline: given an initial parameter set estimate THETA_0, choose V and
% compute K, Y, H_c, H_1_hat, ..., H_n_alpha_hat, H_Q, and H_R

% Choose V matrix

% creates a V matrix with dimensions 2^5 x 2
V = generate_polytope3(2, 5);
% Take first 8 rows of V
V = V(1:8,:);
% append row of F that corresponds to row of 0s in G
V = [V; 0 -3.33];

% Calculate H_c and K, and lambda using lemma7():
[H_c, K, lambda] = lemma7(A0, A1, A2, A3, B0, B1, B2, B3, PI_theta, pi_t, F, G, V);
% display(value(lambda))

% Calculate H_1_hat, ..., H_n_alpha_hat by using lemma8():
H_hat = lemma8(PI_theta, pi_t, A0, A1, A2, A3, B0, B1, B2, B3, K, V);
% Access elements of H_hat by using H_hat{i}


% Calculate H_Q and H_R:

















%% Online Section of the Proposed Algorithm 


% Do the parameter set update with this function to get pi_t_plus_one
pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,u_t_1,x_t,PI_theta,PI_w,pi_t,pi_w);

% calculate vertices of the newly updated parameter set 
vertices = con2vert(PI_theta,(pi_t_plus_one)');
% disp(V)



