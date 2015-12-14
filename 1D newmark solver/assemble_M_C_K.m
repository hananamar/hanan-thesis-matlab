function [M,C,K] = assemble_M_C_K
    % assebles matrices M,C,K (1D)
    global nodes connections ID Nel Nun Nnp BC_S c rho T;

M = zeros(Nun,Nun);
C = zeros(Nun,Nun);
K = zeros(Nun,Nun);

catM = catalogueM;
catK = catalogueK*c^2;

for e = 1:Nel
    e_nodes = connections(:,e);
    
    %x = [nodes(1,e_nodes(1)); nodes(1,e_nodes(2))];
    %y = [nodes(2,e_nodes(1)); nodes(2,e_nodes(2))];
    for i1 = 1:length(e_nodes)
        un_a = ID(connections(i1,e));
        if un_a ~= 0
            for i2 = 1:length(e_nodes)
                un_b = ID(connections(i2,e));
                if un_b ~= 0
                    K(un_a,un_b) = K(un_a,un_b) + catK(i1,i2);
                    M(un_a,un_b) = M(un_a,un_b) + catM(i1,i2);
                    if i1 == i2 && ismember(connections(i1,e),BC_S)
                        if connections(i1,e) == 1
                            C(un_a,un_b) = C(un_a,un_b) + c;
                        end
                        if connections(i1,e) == Nnp
                            C(un_a,un_b) = C(un_a,un_b) + c;
                        end
                    end
                end
            end
        end
    end
end