function pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,u_t_1,x_t,PI_theta,PI_w,pi_t,pi_w)

D = [(A1 * x_t_1+ B1* u_t_1) (A2 * x_t_1+ B2* u_t_1) (A3 * x_t_1+ B3* u_t_1)];
P_t = -PI_w * D;

d_t = A0 * x_t_1 + B0 * u_t_1 - x_t;

Q_t = pi_w + PI_w *d_t;

% do the optimization
% have to do this optimization in a forloop 2p times to get vector
% pi_t_plus_one of length 2p x 1

% set YALMIP decision variables
H_i = sdpvar(1,10);
var_pi = sdpvar(1);


% set options for solver
% options = sdpsettings('solver','gurobi');
options = sdpsettings('solver','gurobi','verbose',0);

for i=1:length(pi_t)
    % set constraints
    % Constraints = [H_i * [PI_theta; P_t] >= PI_theta(1,:)];
    % Constraints = [Constraints, H_i * [PI_theta; P_t] <= PI_theta(1,:)];
    % Constraints = [Constraints, H_i >= 0, H_i * [pi_t; Q_t] <= var_pi];
    % Constraints = [Constraints, var_pi >= 0];

    Constraints = [H_i * [PI_theta; P_t] == PI_theta(i,:), H_i * [pi_t; Q_t] <= var_pi, H_i >= 0];


    % set objective that we want to minimize 
    Objective = var_pi;

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