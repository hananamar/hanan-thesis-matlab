function return_F = F_conc(t)
    global t_0 F F_node ID Nun;
    return_F = zeros(Nun,1);
    if t >= 0 && t < t_0
        return_F(ID(F_node)) = F;
    end
end