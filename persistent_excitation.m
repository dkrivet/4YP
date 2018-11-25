function D = persistent_excitation(predicted_state, predicted_input, A1, A2, A3, B1, B2, B3)

p = 10;

j = 1;
for i = 1:2:2p-2
    D(i:i+1,:) = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, predicted_state(j,:)', predicted_input(j));
    j = j + 1;
end 


