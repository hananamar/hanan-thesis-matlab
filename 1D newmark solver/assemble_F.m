function F = assemble_F(t)
    global Nun BC_N BC_D nodes ID connections bc_q bc_u Nel dt T;
    %assuming boundary condition's time derivatives are 0
    %F_conc_size = 0; %[N]
    f_size = 0;      %[N/m]
    F = zeros(Nun,1);
    
    catM = catalogueM;
    catK = catalogueK;

    for e = 1:Nel
        e_nodes = connections(:,e);
        for i1 = 1:length(e_nodes)
            un_a =  ID(connections(i1,e));
            if un_a ~= 0
                for i2 = 1:length(e_nodes)
                    un_b = ID(connections(i2,e));
                    F(un_a) = F(un_a) + catM(i1,i2)*f_size;
                    if ismember(find(ID == un_a), BC_N)
                        F(un_a) = F(un_a) + catM(i1,i2)*bc_q;
                    end
                    if un_b == 0
                        F(un_a) = F(un_a) - T*catK(i1,i2)*bc_u(round(t/dt) + 1);
                    end
                end
            end
        end
    end
    
    %F = F+F_conc(t);
end