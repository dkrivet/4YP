function D = compute_D_of_x_and_u(A1, A2, A3, B1, B2, B3, x, u)

D = [(A1 * x+ B1* u) (A2 * x+ B2* u) (A3 * x+ B3* u)];

end 