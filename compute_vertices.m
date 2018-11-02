function vertices = compute_vertices(PI_theta, pi_t)


% create the diagonal matrix:
diag_matrix = [1 1 1 1 -1 -1 -1 -1;1 1 -1 -1 1 1 -1 -1;1 -1 1 -1 1 -1 1 -1];

a = length(pi_t);

theta_0 = (1/2) * (pi_t(1:a/2) - pi_t((a/2) + 1: end));

if pi_t(1:a/2) == pi_t((a/2) + 1: end)
    for i = 1:length(diag_matrix(1,:))
        vertices(i,:) = diag_matrix(:,i) .* pi_t(1:a/2);
    end
else
    pi_t = [(1/2)*(pi_t(1:a/2) + pi_t((a/2)+1:end)) ; (1/2)*(pi_t(1:a/2) + pi_t((a/2)+1:end))];
    for i = 1:length(diag_matrix(1,:))
        vertices(i,:) = diag_matrix(:,i) .* pi_t(1:a/2);
        %size(vertices(i,:))
        %size(theta_0)
        vertices(i,:) = vertices(i,:) + theta_0';
    end
end

