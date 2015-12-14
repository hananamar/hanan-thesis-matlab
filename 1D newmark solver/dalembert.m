function analytic = dalembert(time_steps, nodes, u_0)
%calculates dalembert's formula, assuming derivative is 0
global c;
time_steps_count=length(time_steps);
nodes_count=length(nodes);
left_boundary = nodes(1,1);
right_boundary = nodes(1,end);

analytic = zeros(nodes_count, time_steps_count);
    for i=1:time_steps_count
        t = time_steps(i);
        
        for j=1:nodes_count
            node_x = nodes(1, j);
            
            x1 = node_x + c*t;
            x2 = node_x - c*t;
            u_0_interp = interp1(nodes(1,:), u_0', [x1; x2]);
            if x1 < left_boundary || x1 > right_boundary
                u_0_interp(1) = 0;
            end
            if x2 < left_boundary || x2 > right_boundary
                u_0_interp(2) = 0;
            end
            
            analytic(j,i) = 0.5*(u_0_interp(1) + u_0_interp(2));
            if isnan(analytic(j,i))
                bug ;
            end
        end
    end
end