function lambda_t = update_lambda_t(vertices, H_hat)

theta_hat_transpose = [ones(length(vertices(:,1)),1) vertices];

maximum = 0;
for j = 1:length(vertices(:,1))
    for i = 1:length(H_hat)
        temp = theta_hat_transpose(j,:) * H_hat{i} * ones(9,1);
        if temp >= maximum
            maximum = temp;
        end 
    end
end 


lambda_t = maximum;

end 
