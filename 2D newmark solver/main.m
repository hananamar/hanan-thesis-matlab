%% Main routine
% Main routine is used to set parameters, run calculation and draw figures
clc;
%clear all;
close all;

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
global F h betta gamma rho T dt nodes connections ID Nnp t_0 number_of_steps D_node F_node F_position D_position;

%% Define Newmark Method Params
betta = 1/12;
gamma = 1/2;

%% Choose Mesh Parameter
h = 1/20; % size of a single element [m]

%% Choose Time Step Size
dt = 0.02; % length of a single step in time [s]
t_0 = 1;
t_end = 10; % final time [s]
time_steps = 0 : dt : t_end;
number_of_steps = length(time_steps);

%% Physical Parameters
L = 1; %[m]
width = 2*L; %[m]
height = L; %[m]
F = 0;    % [N]
rho = 5; % [Kg/m^2]
T = 4;    % [N/m]

%% Configure Mesh
[nodes,connections,ID] = meshGenerator(h, width, height);

%% Initial Conditions
u_0 = 0;
v_0 = 0;

%% Solve for u
tt = cputime; % start time measurment

solution = NewmarkSolver(u_0, v_0);

disp(['Elapsed time: ' num2str(cputime - tt) ' sec']);% stop time measurment

%% Create Graphs

print_3d_surface_at_times = 1;
if print_3d_surface_at_times
    my_print_3d(solution(1:Nnp, t_0/dt), nodes, connections, 'Displacement at t=t_0');
    my_print_3d(solution(1:Nnp, (t_end/2)/dt), nodes, connections, 'Displacement at t=t_e_n_d/2');
    my_print_3d(solution(1:Nnp, number_of_steps - 1), nodes, connections, 'Displacement at t=t_e_n_d');
end

print_detector_reading = 1;
if print_detector_reading
    my_print_node(detector_readings, time_steps);
end
