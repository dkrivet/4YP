function mu = compute_mu(A1, A2, A3, B1, B2, B3)

gray_code = [1 1 1;1 1 2; 1 2 1; 1 2 2; 2 1 1;2 1 2; 2 2 1; 2 2 2];

x1 = [-100 100];
x2 = [-0.3 100];
u = [-1 1];

for i =1:length(gray_code(:,1))
    vertices(i,:) = [x1(gray_code(i,1)) x2(gray_code(i,2)) u(gray_code(i,3))];
end

% each row of vertices containt [x1_i x2_i u_i]

lambda_bar = 0;

for i =1:length(vertices(:,1))
    x = [vertices(i,1); vertices(i,2)];
    u = vertices(i,3);
    D = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x, u);
    temp = max(eig(D' * D));
    if temp > lambda_bar
        lambda_bar = temp;
    end
end

mu = 1/lambda_bar;

end
    

