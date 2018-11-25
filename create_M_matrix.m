function M = create_M_matrix(N, element)

M = ones(2 * N,2);
for i = 1:2:2 * N
    M(i:i+1,:) = element ^ (ceil(i/2));
end