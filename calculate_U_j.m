function U_j = calculate_U_j(R_j, V)

% for i = 1:length(R_j)
%     U_j{i} = [];
% end

U_j = zeros(2, length(V(:,1)), length(V(:,1)));

I = eye(length(V(:,1)));

for i = 1:length(R_j)
    U_j(:,:,i) = (V([R_j(i,1), R_j(i,2)],:))^(-1) * I([R_j(i,1), R_j(i,2)],:);
end 

