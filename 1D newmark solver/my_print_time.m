function my_print_time(v, x_values, time)
    % prints a displacement at all nodes in a certain time
    % v is vector of values at all nodes at given time step
    
    figure;
    hold on;
    grid on;   
    
    title(['Displacement at Time = ', num2str(time), '[sec]']);
    xlabel('X [m]');
    ylabel(['Displacement [m]']);
    
    plot(x_values, v);
end