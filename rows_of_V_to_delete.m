function indeces = rows_of_V_to_delete(alpha_bar)

indeces = [];

for i = 1:length(alpha_bar)
    if alpha_bar(i) <= 1
        indeces = [indeces, i];
    end
end

end
