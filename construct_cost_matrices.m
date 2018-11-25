function [H, f_transpose] = construct_cost_matrices(A0, A1, A2, A3, B0, B1, B2, B3, K, N, R, Q, x_k, theta_hat)

% The equation defining PHI should have A_hat and B_hat, not A0 and B0
% ask Prof Cannon what A_hat and B_hat is 
% A_hat = A0;
% B_hat = B0;

[A_hat, B_hat] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,theta_hat);
PHI = A_hat + B_hat * K;



% Create M matrix
% M = ones(2 * N,2);
% for i = 1:2:2 * N
%     M(i:i+1,:) = PHI ^ (ceil(i/2));
% end
M = create_M_matrix(N, PHI);

% Create C matrix
a = zeros(2,1);
% 
% % This surely can be done in a better way:
% C(:,1) = [B_hat ; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat; PHI^4 * B_hat; PHI^5 * B_hat; PHI^6 * B_hat; PHI^7 * B_hat; PHI^8 * B_hat; PHI^9 * B_hat];
% C(:,2) = [a ; B_hat; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat; PHI^4 * B_hat; PHI^5 * B_hat; PHI^6 * B_hat; PHI^7 * B_hat; PHI^8 * B_hat];
% C(:,3) = [a ; a; B_hat; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat; PHI^4 * B_hat; PHI^5 * B_hat; PHI^6 * B_hat; PHI^7 * B_hat];
% C(:,4) = [a ; a; a; B_hat; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat; PHI^4 * B_hat; PHI^5 * B_hat; PHI^6 * B_hat];
% C(:,5) = [a ; a; a; a; B_hat; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat; PHI^4 * B_hat; PHI^5 * B_hat];
% C(:,6) = [a ; a; a; a; a; B_hat; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat; PHI^4 * B_hat];
% C(:,7) = [a ; a; a; a; a; a; B_hat; PHI * B_hat; PHI^2 * B_hat; PHI^3 * B_hat];
% C(:,8) = [a ; a; a; a; a; a; a; B_hat; PHI * B_hat; PHI^2 * B_hat];
% C(:,9) = [a ; a; a; a; a; a; a; a; B_hat; PHI * B_hat];
% C(:,10) = [a ; a; a; a; a; a; a; a; a; B_hat];

C = zeros(2 * N, N);
first_row = [B_hat a a a a a a a a a];
C = create_C_matrix(C, first_row, 1, N, PHI, B_hat, 0);




% Create R_bar
R_bar = zeros(10,10);

for i = 1:10
    for j = 1:10
        if i == j
            R_bar(i,j) = R;
        end
    end
end

% Calculate P matrix: solution to the equation: 
% P - (A_hat + B_hat * K)' * P * (A_hat + B_hat * K) = Q + K' * R * K
P = dlyap((A_hat + B_hat * K)', -Q - K' * R * K);
P = -P;

% Create Q_bar
Q_bar = zeros(2 * N, 2 * N);
for i = 1:2:2*N
    for j = 1:2:2*N
        if (i == j) && i ~= 2 * N -1
            Q_bar(i:i+1,j:j+1) = Q;
        elseif (i == j) && i == 2 * N - 1
            Q_bar(i:i+1,j:j+1) = P;
        end
    end
end 

            




% values to return:
H = C' * Q_bar * C + R_bar;
f_transpose = x_k' * M' * Q_bar * C;

end