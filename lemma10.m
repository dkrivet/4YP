function [H_Q, H_R] = lemma10(V, Q, R, K)

% find H_Q:

% set YALMIP decision variables
H = sdpvar(1,length(V(:,1)));

% set options for solver
options = sdpsettings('solver','gurobi');

% set objective
Objective = H * ones(length(V(:,1)), 1);

for i = 1:length(Q(:,1))
    % set constraints:
    Constraints = [H >= 0, H * V == Q(i,:)];

    % solve the optimization 
    sol = optimize(Constraints, Objective, options);
    
    H_Q(i) = value(H) * (ones(length(V(:,1)), 1));
end



% find H_R:
Constraints = [H >= 0, H * V == R * K];

% solve the optimization 
sol = optimize(Constraints, Objective, options);

H_R = H * ones(length(V(:,1)), 1);

