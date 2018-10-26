function pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,u_t_1,x_t,PI_theta,PI_w,pi_t,pi_w)
% A0 = [0.5 0.2;-0.1 0.6];
% A1 = [0.042 0;0.072 0.03];
% A2 = [0.015 0.019;0.009 0.035];
% A3 = [0 0;0 0];
% B0 = [0; 0.5];
% B1 = [0;0];
% B2 = [0;0];
% B3 = [0.0397;0.059];

% Define the initial condition:
% x_t_1 = [3; 6];

% These following two values I'm not sure how to obtain so have chosen
% random values
% u_t_1 = -1;  % define arbitrary initial control input

% what is the next state
% x_t = [2.9; 3.2];


% PI_theta = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
% PI_w = [1 0;0 1;-1 0;0 -1];
D = [(A1 * x_t_1+ B1* u_t_1) (A2 * x_t_1+ B2* u_t_1) (A3 * x_t_1+ B3* u_t_1)];
P_t = -PI_w * D;


% pi_t = [1; 1; 1; 1; 1; 1];
% pi_w = [0.1; 0.1; 0.1; 0.1];
d_t = A0 * x_t_1 + B0 * u_t_1 - x_t;

Q_t = pi_w + PI_w *d_t;

% do the optimization
% have to do this optimization in a forloop 2p times to get vector
% pi_t_plus_one of length 2p x 1

% set YALMIP decision variables
H_i = sdpvar(1,10);
var_pi = sdpvar(1);

for i=1:length(pi_t)
    % set constraints
    % Constraints = [H_i * [PI_theta; P_t] >= PI_theta(1,:)];
    % Constraints = [Constraints, H_i * [PI_theta; P_t] <= PI_theta(1,:)];
    % Constraints = [Constraints, H_i >= 0, H_i * [pi_t; Q_t] <= var_pi];
    % Constraints = [Constraints, var_pi >= 0];

    Constraints = [H_i * [PI_theta; P_t] == PI_theta(i,:), H_i * [pi_t; Q_t] <= var_pi, H_i >= 0];


    % set objective that we want to minimize 
    Objective = var_pi;

    % set options for solver
    % options = sdpsettings('verbose',1,'solver','gurobi');
    options = sdpsettings('solver','gurobi');

    % Solve the optimization
    % sol = optimize([Constraints, -1000<= var_pi <= 1000], Objective);
    sol = optimize(Constraints, Objective, options);

    % This is the value of pi_t_plus_one
    
    % disp(value(var_pi))
    pi_t_plus_one(i) = value(var_pi);
    
end

% disp(pi_t_plus_one)
% This code is to print the convex polytopic set
% theta = sdpvar(3,1);
% plot(PI_theta * theta <= pi_t)
% con2vert(PI_theta,pi_t)



% after computing pi_t_plus_one plot it:
% plot(PI_theta * theta <= (pi_t_plus_one)')


% computes vertices of polytope given H-form polytope
% con2vert(PI_theta,(pi_t_plus_one)')
end