function [predicted_state, predicted_input] = compute_predicted_state_and_input(predicted_v, K, x_t, N, point_estimate, A0, A1, A2, A3, B0, B1, B2, B3)

predicted_state = zeros(N,1);
for i = 1:N
    predicted_state(i) = K * x_t + predicted_v(i);
end 

[A, B] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,point_estimate);

M = zeros(2 * N, 2)

for i = 1:2:2*N
    M(i:i+1,:) = A^(ceil(i/2));
end




    

