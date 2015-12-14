function [M,K] = assemble_M_K
    % assebles matrices M,K
    global nodes connections ID T rho;
Nel = length(connections);
Nun = max(ID);

M = zeros(Nun,Nun);
K = zeros(Nun,Nun);

for e = 1:Nel
    e_nodes = connections(:,e);
    
    x = [nodes(1,e_nodes(1)); nodes(1,e_nodes(2)); nodes(1,e_nodes(3))];
    y = [nodes(2,e_nodes(1)); nodes(2,e_nodes(2)); nodes(2,e_nodes(3))];
    for i1 = 1:length(e_nodes)
        un_a = ID(connections(i1,e));
        if un_a ~= 0
            for i2 = 1:length(e_nodes)
                un_b = ID(connections(i2,e));
                catM = catalogueM(x,y);
                if un_b ~= 0
                    K(un_a,un_b) = K(un_a,un_b) + T*catalogueK(i1,i2,x,y);
                    M(un_a,un_b) = M(un_a,un_b) + rho*catM(i1,i2);
                end
            end
        end
    end
end