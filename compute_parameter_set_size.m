function radial_size = compute_parameter_set_size(pi_t)

pi_1 = pi_t(1:length(pi_t)/2);
pi_2 = -pi_t(length(pi_t)/2 + 1: end);

radial_size = 0;

for i = 1: length(pi_1)
    temp_distance = abs(pi_1(i) - pi_2(i));
    if temp_distance > radial_size
        radial_size = temp_distance;
    end
end


end