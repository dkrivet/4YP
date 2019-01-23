function [feasible_set_x, feasible_set_y, infeasible_set_x, infeasible_set_y] = compute_feasible_set()

x_range = linspace(-10, 10, 25);
y_range = linspace(-0.3, 22, 40);
feasible_set_x = [];
feasible_set_y = [];
infeasible_set_x = [];
infeasible_set_y = [];


for i = 1:length(x_range)
    for j = 1:length(y_range)
        sol = main([x_range(i);y_range(j)]);
        if sol.problem == 0
            feasible_set_x = [feasible_set_x, x_range(i)];
            feasible_set_y = [feasible_set_y, y_range(j)];
        else
            infeasible_set_x = [infeasible_set_x, x_range(i)];
            infeasible_set_y = [infeasible_set_y, y_range(j)];
        end
        i
        j
    end
end


scatter(feasible_set_x, feasible_set_y, 25, 'b')
hold on
scatter(infeasible_set_x, infeasible_set_y, 25, 'r')

end

