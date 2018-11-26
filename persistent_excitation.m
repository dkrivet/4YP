function du_to_return = persistent_excitation(predicted_state, predicted_input, A1, A2, A3, B1, B2, B3, var_beta)

p = 10;

j = 1;
for i = 1:2:2*p
    D(i:i+1,:) = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, predicted_state(j,:)', predicted_input(j));
    j = j + 1;
end


du = sdpvar(10,1);

% have to define du (I think as sdpvar)
j = 1;
for i = 1:2:2*p
    L_of_du(i:i+1,:) = compute_L_of_du(B1, B2, B3, du(j,1));
    j = j + 1;
end 

options = sdpsettings('solver','gurobi','verbose',1);

Objective = du;

I = eye(3);
Constraints = [D' * D + L_of_du' * D + D' * L_of_du >= p * var_beta^2 * I];

% solve the optimization
sol = optimize(Constraints, Objective, options);

du_to_return = value(du);


end 


