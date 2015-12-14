function norm = calc_norm(v)
    sum = 0;
    for j = 1:length(v)
        sum = sum + v(j)^2;
    end
    norm = sqrt(sum/length(v));
end