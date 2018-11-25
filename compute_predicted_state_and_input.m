function [predicted_state, predicted_input] = compute_predicted_state_and_input(predicted_v, K, x_t, N, point_estimate, A0, A1, A2, A3, B0, B1, B2, B3)

predicted_input = zeros(N,1);
for i = 1:N
    predicted_input(i) = K * x_t + predicted_v(i);
end 

[A, B] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,point_estimate);

M = create_M_matrix(N, A);


a = zeros(2,1);
C = zeros(2 * N, N);
first_row = [B a a a a a a a a a];
C = create_C_matrix(C, first_row, 1, N, A, B, 0);


% intiailize predicted state so that row i holds value x_i_given_k
predicted_state = zeros(N, 2);
for i = 1:2:2 * N
    predicted_state(ceil(i/2),:) = (M(i:i+1,:) * x_t + C(i:i+1,:) * predicted_input)';
end 


    

