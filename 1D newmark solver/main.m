%% Main routine
% Main routine is used to set parameters, run calculation and draw figures
clc;
clear all;
%close all;

%% Global parameters:
% F = concentrated force value
% F_node = concentrated force node number
% t_0 = duration of concentrated force
% h = mesh parameter
% Nnp = total number of nodes
% betta = Newmark method parameter
% gamma = Newmark method parameter
% rho = physical attribute of material rho
% T = physical attribute of material T
% t_end = final time
% nodes = mesh array of nodes and their respective coordinates
% connections = mesh array of elements and their respective nodes
% ID = mesh array of nodes and their corresponding unknown variables
% dt = length of a single step in time
global h rho T c dt nodes connections ID Nnp t_0 number_of_steps BC_D BC_N BC_S bc_q bc_u alpha gamma betta;

%% Define HHT Method Parameter
alpha = 0;

%% Define Newmark Method Params
gamma = (1-2*alpha)/2;
betta = (1-alpha)^2/4;

%% Physical Parameters
L = 1; %[m]
width = L; %[m]
rho = 5; % [Kg/m^2]
T = 5;    % [N/m]
c = sqrt(T/rho);

%% Choose Time Step Size
dt = 0.02; % length of a single step in time [s]
t_0 = 0;
number_of_steps = 51;
t_end = dt*(number_of_steps-1); % final time [s]
time_steps = 0 : dt : t_end;

%% Boundary Conditions
% array of boundary conditions type, like so: [left right]
% 1 = closed (dirichlet), 0 = open (neumann), 2 = no reflection (sommerfeld)
boundary_conditions = [2 2];
bc_u = zeros(1,number_of_steps);
% m=1;
% first=1;
% for i=1:number_of_steps
%     if sin(time_steps(i)*pi) >= 0 && first
%         m=i;
%     else
%         first = 0;
%     end
% end
% bc_u(1:m) = sin(time_steps(1:m)*pi);
bc_q = 0;

%% Configure Mesh
h = 1/30; % size of a single element [m]
[nodes, connections, ID, BC_D, BC_N, BC_S] = meshGenerator1D(h, width, boundary_conditions);

%% Initial Conditions - sinusiodal
x = linspace(-pi/2,3/2*pi,ceil(Nnp/3));
y = sin(x)+1;
% dy = cos(x+pi/2);

u_0 = zeros(Nnp,1);
v_0 = zeros(Nnp,1);
u_0(ceil(Nnp/3):2*ceil(Nnp/3)-1) = y;
% v_0(ceil(Nnp/3):2*ceil(Nnp/3)-1) = dy;

for i=1:length(BC_D)
    if BC_D(i) > 0
        u_0(BC_D(i)) = bc_u(i);
    end
end

%% Calculate analytical solution
analytic = dalembert(time_steps, nodes, u_0); 

%% Solve for u
tt = cputime; % start time measurment

%solution = NewmarkSolver(u_0, v_0);
%% insert func% Newmark solver
% required params:
    % argument assemble_M_K_F should be a function handle for assembling matrices M,K
    % argument u_0 should be the initial condition of the solution
    % argument v_0 should be the initial condition of the first time derivative of the solution
    % argument u_boundary is the boundary condition
% returns a matrix of Nnp x number_of_steps with solution
    solution = zeros(Nnp * 2, number_of_steps); % column i of array "solution" is values of solution in time step i-1
    % first half of rows in 'solution' are displacement at open nodes
    % second half of rows in 'solution' are velocity at open nodes
    
    % setting initial condition values (t=0)
    solution(1:Nnp, 1) = u_0;
    solution(Nnp+1:2*Nnp, 1) = v_0;
    
    % setting uniform boundary conditions
%     for i=1:length(BC_D)
%         solution(BC_D(i), :) = bc_u;
%     end
    % setting changing boundary conditions
    if ~isempty(BC_D)
        for i = 1:length(BC_D)
            solution(BC_D(i), :) = bc_u;
        end
    end
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
    
    [M,C,K] = assemble_M_C_K(); % assemble matrices
    M_star = M + (1+alpha) * gamma * dt * C + (1+alpha) * betta * dt^2 * K; 
        
    % initialize a_n at t=0
    F_0 = assemble_F(0);
    a_n = M_star\(F_0 - C * variables_v_n - K * variables_d_n);
    
    %% iterations phase
    F_t_n = F_0;
    F_t_n_plus_1 = 0;
    for n = 1 : number_of_steps-1
        % solve Newmark
        % t_n = dt * (n-1);
        t_n_plus_1 = dt * n;
        d_n_plus_1_predicted = variables_d_n + dt * variables_v_n + (1-2*betta)/2 * dt^2 * a_n;
        v_n_plus_1_predicted = variables_v_n + dt * (1-gamma) * a_n;
        F_t_n_plus_1 = assemble_F(t_n_plus_1);
        F_star = (1+alpha)*F_t_n_plus_1 - alpha*F_t_n - C*( (1+alpha)*v_n_plus_1_predicted - alpha*variables_v_n ) - K*( (1+alpha)*d_n_plus_1_predicted - alpha*variables_d_n );
        
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
        F_t_n = F_t_n_plus_1;
    end
    
%% insert func end

disp(['Elapsed time: ' num2str(cputime - tt) ' sec']);% stop time measurment

%% Create Graphs

print_3d_surface_at_times = 1;
if print_3d_surface_at_times
    my_print_2d_and_time(solution(1:Nnp, :), time_steps, nodes, connections, 'Displacement');
    %my_print_2d_and_time(analytic, time_steps, nodes, connections, 'Displacement');
end

print_error = 1;
if print_error
    errors = zeros(1,number_of_steps);
    for i=1:number_of_steps
        for j = 1:Nnp
            errors(i) = errors(i) + (analytic(j,i) - solution(j,i))^2;
        end
        errors(i) = sqrt(errors(i));
    end
    
    my_print_node(errors, 0:dt:t_end);
end
