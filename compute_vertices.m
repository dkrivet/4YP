function vertices = compute_vertices(PI_theta, pi_t)
rows_of_PI_theta = size(PI_theta);
rows_of_PI_theta = rows_of_PI_theta(1);

% create the diagonal matrix:
diag_matrix = ones(rows_of_PI_theta, rows_of_PI_theta);
diag_matrix(:,end) = -1;

for i = 2:rows_of_PI_theta
    diag_matrix(i,i) = -1;
end


A = diag_matrix * pi_t;


% delete non-unique rows of A
B = unique(A,'rows')