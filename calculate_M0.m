function M0 = calculate_M0(x, u, i, A1, A2, A3, B1, B2, B3)

M0 = 0;

if i < 10
    for j = 1:i
        D_k = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x(:,j), u(j));
        M0 = M0 + D_k' * D_k; 
    end
else
    for j = i-9:i
        D_k = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x(:,j), u(j));
        M0 = M0 + D_k' * D_k;        
    end
end



end