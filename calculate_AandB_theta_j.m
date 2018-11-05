function [A_theta_j, B_theta_j] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,vertex)

A_theta_j = A0 + A1 * vertex(1) + A2 * vertex(2) + A3 * vertex(3);
B_theta_j = B0 + B1 * vertex(1) + B2 * vertex(2) + B3 * vertex(3);

end 