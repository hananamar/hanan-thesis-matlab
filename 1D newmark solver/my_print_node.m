function my_print_node(v, time_steps)
    % prints a time vs. temperature graph at origin
    % v is vector of values at node for each time step
    
    figure;
    hold on;
    grid on;   
    
    title('Displacement error through time');
    xlabel('Time [sec]');
    ylabel(['Displacement error size [m]']);
    
    plot(time_steps, v);
end