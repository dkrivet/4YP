function w_bar = compute_w_bar(PI_w, pi_w, V)

% set YALMIP decision variables
w = sdpvar(length(PI_w(1,:)), 1);

% set constraints
Constraints = [PI_w * w <= pi_w];

% set options for solver
options = sdpsettings('solver','gurobi','verbose',0);

for i = 1:length(V(:,1))
    % set objective
    Objective = - V(i,:) * w;
    
    % solve the optimization
    sol = optimize(Constraints, Objective, options);
    
    w_bar(i) = value(V(i,:) * w);
end
end 
