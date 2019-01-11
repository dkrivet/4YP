function R_j = compute_R_j(V, tube_vertices, alpha)

no_of_vertices = length(tube_vertices(:,1));

for i = 1:no_of_vertices
    R_j{i} = [];
end 


for j = 1:no_of_vertices
    for r = 1:length(tube_vertices(:,1))
        if ((alpha(r) - 0.0001) < V(r,:) * tube_vertices(j,:)') && (alpha(r) + 0.0001 > V(r,:) * tube_vertices(j,:)')
            R_j{j} = [R_j{j}, r];
        end
    end
end 

end