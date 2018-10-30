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


%% Offline Section of the Proposed Algorithm 
% Offline: given an initial parameter set estimate THETA_0, choose V and
% compute K, and lambda


















%% Online Section of the Proposed Algorithm 


% Do the parameter set update with this function to get pi_t_plus_one
pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,u_t_1,x_t,PI_theta,PI_w,pi_t,pi_w);

% calculate vertices of the newly updated parameter set 
[V,nr] = con2vert(PI_theta,(pi_t_plus_one)');
% disp(V)



