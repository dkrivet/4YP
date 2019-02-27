function pi_t_plus_one = varying_parameter_set_update(A0,A1,A2,A3,B0,B1,B2,B3,x_t_1,u_t_1,x_t,PI_theta,PI_w,pi_t,pi_w,U,h)

D = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x_t_1, u_t_1);
d = A0 * x_t_1 + B0 * u_t_1 - x_t;


theta_star = sdpvar(3,1);
theta_tilda = sdpvar(3,1);
w = sdpvar(2,1);

Constraints = [PI_theta * theta_star <= pi_t];
Constraints = [Constraints, -d - D * theta_star == D * theta_tilda + w];
Constraints = [Constraints, U * theta_tilda <= h];
Constraints = [Constraints, PI_w * w <= pi_w];

options = sdpsettings('solver','gurobi','verbose',1);

for i = 1:length(PI_theta(:,1))
    Objective = -PI_theta(i,:) * theta_star;
    
    sol = optimize(Constraints, Objective, options);
    
    if sol.problem ~= 0
        fprintf('ATTENTION!! \n Varying Parameter Set Update Optimization had a problem \n')
    end
    
    pi_t_plus_one(i) = value(PI_theta(i,:) * theta_star);

end