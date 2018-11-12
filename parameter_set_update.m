function pi_t_plus_one = parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,u_t_1,x_t,PI_theta,PI_w,pi_t,pi_w)

% D = [(A1 * x_t_1+ B1* u_t_1) (A2 * x_t_1+ B2* u_t_1) (A3 * x_t_1+ B3* u_t_1)];
D = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_t_1, u_t_1);
P_t = -PI_w * D;

d_t = A0 * x_t_1 + B0 * u_t_1 - x_t;

Q_t = pi_w + PI_w *d_t;

% do the optimization
% have to do this optimization in a forloop 2p times to get vector
% pi_t_plus_one of length 2p x 1

% set YALMIP decision variables
H_i_1 = sdpvar(1,10);
H_i_2 = sdpvar(1,10);
var_pi_1 = sdpvar(1);
var_pi_2 = sdpvar(1);


% set options for solver
% options = sdpsettings('solver','gurobi');
options = sdpsettings('solver','gurobi','verbose',0);

for i=1:length(pi_t)/2
    % set constraints
    Constraints = [H_i_1 * [PI_theta; P_t] == PI_theta(i,:), 
        H_i_2 * [PI_theta; P_t] == PI_theta(i+3,:), 
        H_i_1 * [pi_t; Q_t] <= var_pi_1, 
        H_i_2 * [pi_t; Q_t] <= var_pi_2,
        H_i_1 >= 0, 
        H_i_2 >= 0];


    % set objective that we want to minimize 
    Objective = var_pi_1 + var_pi_2;

    % Solve the optimization
    % sol = optimize([Constraints, -1000<= var_pi_1 <= 1000, -1000<= var_pi_2 <= 1000], Objective);
    sol = optimize(Constraints, Objective, options);

    % This is the value of pi_t_plus_one
    pi_t_plus_one(i) = value(var_pi_1);
    pi_t_plus_one(i+3) = value(var_pi_2);
    
end

end