function mu = compute_mu(A1, A2, A3, B1, B2, B3, F, G)

% set YALMIP decision variables 
x = sdpvar(2,1);
u = sdpvar(1);

% options = sdpsettings('solver','gurobi');
% options = sdpsettings('solver','gurobi','verbose',0);

% set objective
D = [(A1 * x + B1 * u) (A2 * x + B2 * u) (A3 * x + B3 * u)];
% test = @(x,u) [(A1*x+B1*u)'*(A1*x+B1*u) (A1*x+B1*u)'*(A2*x+B2*u) (A1*x+B1*u)'*(A3*x+B3*u);(A2*x+B2*u)'*(A1*x+B1*u) (A2*x+B2*u)'*(A2*x+B2*u) (A2*x+B2*u)'*(A3*x+B3*u);(A3*x+B3*u)'*(A1*x+B1*u) (A3*x+B3*u)'*(A2*x+B2*u) (A3*x+B3*u)'*(A3*x+B3*u)];
% max_eig_value = max(eig(test(x_opt,u_opt)));
% Objective = -max_eig_value; 
Objective = -norm(D);
% set constraints
Constraints = [F * x + G * u <= ones(2,1)];

% solve the optimization
sol = optimize(Constraints, Objective);


mu = value(norm(D));
% mu = 1/(-max(eig(D' * D)));

