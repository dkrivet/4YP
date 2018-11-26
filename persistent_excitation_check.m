function bool_matrix = persistent_excitation_check(predicted_state, predicted_input, A1, A2, A3, B1, B2, B3, var_beta)

p = 10;

j = 1;
for i = 1:2:2*p
    D(i:i+1,:) = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, predicted_state(j,:)', predicted_input(j));
    j = j + 1;
end

bool_matrix = (D' * D >= 10 * var_beta^2 * eye(3));


end 


