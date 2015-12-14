function my_print_3d(v, nodes, connections, name)
    global crack_nodes crack_elements;
    % prints a 3d surface of the given solution
    % v is vector of values in each node
    % nodes and connections are arrays from mesh generator
    figure
    hold on
    grid on
    
    title(name)
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('Displacement [m]')
    for i=1:length(nodes)
        XYZ(i,:) = [nodes(1,i), nodes(2,i), v(i)];
    end
    for i=1:length(connections)
        el_x = [XYZ(connections(1,i),1), XYZ(connections(2,i),1), XYZ(connections(3,i),1)];
        el_y = [XYZ(connections(1,i),2), XYZ(connections(2,i),2), XYZ(connections(3,i),2)];
        el_u = [XYZ(connections(1,i),3), XYZ(connections(2,i),3), XYZ(connections(3,i),3)];
        
        color = 'w';
        
        % give crack elements different color
        for j = 1:length(crack_nodes)
            if connections(1,i) == crack_nodes(j,1) || connections(2,i) == crack_nodes(j,1) || connections(3,i) == crack_nodes(j,1)
                color = 'r';
            end
            if connections(1,i) == crack_nodes(j,2) || connections(2,i) == crack_nodes(j,2) || connections(3,i) == crack_nodes(j,2)
                color = 'g';
            end
        end
        fill3(el_x,el_y,el_u,color);
    end
end