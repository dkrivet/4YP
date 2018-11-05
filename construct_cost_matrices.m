function [H, f_transpose] = construct_cost_matrices(A0, A1, A2, A3, B0, B1, B2, B3, K, N, R, Q, x_k, theta_hat)

% The equation defining PHI should have A_hat and B_hat, not A0 and B0
% ask Prof Cannon what A_hat and B_hat is 
% A_hat = A0;
% B_hat = B0;

[A_hat, B_hat] = calculate_AandB_theta_j(B0,B1,B2,B3,A0,A1,A2,A3,theta_hat);
PHI = A_hat + B_hat * K;



% Create M matrix
M = ones(2 * N,2);
for i = 1:2:2 * N
    M(i:i+1,:) = PHI ^ i;
end

% Create C matrix
a = zeros(2,1);
% 
% % This surely can be done in a better way:
% C(:,1) = [B ; PHI * B; PHI^2 * B; PHI^3 * B; PHI^4 * B; PHI^5 * B; PHI^6 * B; PHI^7 * B; PHI^8 * B; PHI^9 * B];
% C(:,2) = [a ; B; PHI * B; PHI^2 * B; PHI^3 * B; PHI^4 * B; PHI^5 * B; PHI^6 * B; PHI^7 * B; PHI^8 * B];
% C(:,3) = [a ; a; B; PHI * B; PHI^2 * B; PHI^3 * B; PHI^4 * B; PHI^5 * B; PHI^6 * B; PHI^7 * B];
% C(:,4) = [a ; a; a; B; PHI * B; PHI^2 * B; PHI^3 * B; PHI^4 * B; PHI^5 * B; PHI^6 * B];
% C(:,5) = [a ; a; a; a; B; PHI * B; PHI^2 * B; PHI^3 * B; PHI^4 * B; PHI^5 * B];
% C(:,6) = [a ; a; a; a; a; B; PHI * B; PHI^2 * B; PHI^3 * B; PHI^4 * B];
% C(:,7) = [a ; a; a; a; a; a; B; PHI * B; PHI^2 * B; PHI^3 * B];
% C(:,8) = [a ; a; a; a; a; a; a; B; PHI * B; PHI^2 * B];
% C(:,9) = [a ; a; a; a; a; a; a; a; B; PHI * B];
% C(:,10) = [a ; a; a; a; a; a; a; a; a; B];
C = zeros(2 * N, N);
for i = 1:2:2*N
    for j = 1:N
        if (j >= i && j ~= 1) || ((abs(j - i) <= floor(j/2)) && j > 3)
            C(i:i+1,j) = a;
        elseif j == ceil(i/2);
            C(i:i+1,j) = B_hat;
        else
            C(i:i+1, j) = A_hat^(ceil(i/2) - j) * B_hat;
        end
    end
end 



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
Q_bar = ones(2 * N, 2 * N);
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


