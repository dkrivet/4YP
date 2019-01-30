function R_j = compute_R_j(V, tube_vertices, alpha)

no_of_vertices = length(tube_vertices(:,1));

% for i = 1:no_of_vertices
%     R_j{i} = [];
% end 
R_j = zeros(no_of_vertices,2);


for j = 1:no_of_vertices
    row_to_append = [];
    for r = 1:length(tube_vertices(:,1))
        if ((alpha(r) - 0.0001) < V(r,:) * tube_vertices(j,:)') && (alpha(r) + 0.0001 > V(r,:) * tube_vertices(j,:)')
            row_to_append = [row_to_append, r];
        end
    end
    R_j(j,:) = row_to_append;
end 

end