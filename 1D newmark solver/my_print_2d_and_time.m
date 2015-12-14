function my_print_2d_and_time(solution, time_vector, nodes, connections, name)
    % prints a 3d surface of the given solution
    % v is vector of values in each node
    % nodes and connections are arrays from mesh generator
    figure
    hold on
    grid on
    
    title(name)
    xlabel('x [m]')
    ylabel('t [sec]')
    zlabel('Displacement [m]')
    for i=1:length(connections)
        for j=1:length(time_vector)-1
            el_x = [nodes(1,connections(1,i)), nodes(1,connections(2,i)), nodes(1,connections(2,i)), nodes(1,connections(1,i))];
            el_y = [time_vector(j), time_vector(j), time_vector(j+1), time_vector(j+1)];
            el_u = [solution(i,j), solution(i+1,j), solution(i+1,j+1), solution(i,j+1)];

            color = 'w';

            fill3(el_x,el_y,el_u,color);
        end
    end
end