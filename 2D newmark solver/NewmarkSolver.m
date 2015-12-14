% Newmark solver
% required params:
    % argument assemble_M_K_F should be a function handle for assembling matrices M,K
    % argument u_0 should be the initial condition of the solution
    % argument v_0 should be the initial condition of the first time derivative of the solution
    % argument u_boundary is the boundary condition
% returns a matrix of Nnp x number_of_steps with solution
function solution = NewmarkSolver(u_0, v_0)
    global Nnp betta gamma dt ID number_of_steps;
    solution = zeros(Nnp * 2, number_of_steps); % column i of array "solution" is values of solution in time step i-1
    % first half of rows in 'solution' are displacement at open nodes
    % second half of rows in 'solution' are velocity at open nodes
    
    % setting initial condition values (t=o)
    solution(1:Nnp, 1) = ones(Nnp, 1) * u_0;
    solution(Nnp+1:2*Nnp, 1) = ones(Nnp, 1) * v_0;
    
    % setting boundary conditions
    % we will skip setting this part due to the fact that boundary conditions are
    % identical to initial conditions.
    
    % initialize array that stores d only in the open nodes
    dof = max(ID);
    variables_d_n = zeros(dof, 1); % displacement variables vector on current step
    variables_v_n = zeros(dof, 1); % velocity variables vector on current step
    for i = 1:Nnp
        if ID(i) > 0
            variables_d_n(ID(i)) = solution(i,1);
            variables_v_n(ID(i)) = solution(i + Nnp,1);
        end
    end
    
    [M,K] = assemble_M_K(); % assemble matrices
    M_star = M + betta * dt^2 * K; % no friction (C = 0)
        
    % initialize a_n at t=0
    a_n = M_star\(F_conc(0) - K * variables_d_n);
    
    %% iterations phase
    for n = 1 : number_of_steps-1
        % solve Newmark
        % t_n = dt * (n-1);
        t_n_plus_1 = dt * n;
        d_n_plus_1_predicted = variables_d_n + dt * variables_v_n + (1-2*betta) * dt^2/2 * a_n;
        v_n_plus_1_predicted = variables_v_n + dt * (1-gamma) * a_n;
        F_star = F_conc(t_n_plus_1) - K * d_n_plus_1_predicted; % no friction (C = 0)
        
        a_n_plus_1 = M_star\F_star;
        variables_d_n_plus_1 = d_n_plus_1_predicted + betta * dt^2 * a_n_plus_1;
        variables_v_n_plus_1 = v_n_plus_1_predicted + gamma * dt * a_n_plus_1;
        
        % store calculated variables in solution array
        for i = 1:dof
            solution(find(ID == i),n+1) = variables_d_n_plus_1(i);
            solution(find(ID == i) + Nnp,n+1) = variables_v_n_plus_1(i);
        end
        
        % next step
        variables_d_n = variables_d_n_plus_1;
        variables_v_n = variables_v_n_plus_1;
        a_n = a_n_plus_1;
    end
    
end