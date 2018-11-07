function C = create_C_matrix(C, first_row, i, N, PHI, B_hat, counter)

C(i:i+1,:) = first_row;

first_row = PHI * first_row;
first_row(1:2,counter+2) = B_hat;
counter = counter + 1;

i = i + 2;

if i < 2 * N
    C = create_C_matrix(C, first_row, i, N, PHI, B_hat, counter);
end



end