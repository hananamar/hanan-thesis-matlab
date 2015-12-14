% Creates mesh properties. Specific for this project's geometry.
% - input 'h' is length of a single element in meters. preferably should be a simple fraction.
% - input 'length' is width of domain in meters.
% - input 'closed_boundaries' represents (boolean) which boundaries are closed: [left, right].
% - output 'nodes' is array of coordinate (x) of each node, and a closed/open code.
% - output 'connections' is array of corresponding nodes for each element. Identical to IEN array.
% - output 'ID' is array of open variables.

function [nodes,connections,ID,BC_D,BC_N,BC_S] = meshGenerator1D(h, width, closed_boundaries)
    global Nnp Nun Nel;
    Nnp = width / h + 1;
    Nel = width / h;
    nodes = zeros(2, Nnp); % initialize
    connections = zeros(2, Nel); % initialize
    ID = zeros(1,Nnp); % initialize
    BC_D = []; % initialize
    BC_N = []; % initialize
    BC_S = []; % initialize

    for i = 1 : Nnp
        nodes(1,i) = (i-1) * h; %populating 'nodes' array
        nodes(2,i) = (closed_boundaries(1)==1 && i == 1) || (closed_boundaries(2)==1 && i == Nnp); % code for open/closed node
    end
    
    for node=1:Nnp-1  %populating 'connections' array
        % 'top' shape elements
        connections(1,node) = node;
        connections(2,node) = node + 1;
    end
    
    counter = 1;
    for i = 1:Nnp
        if nodes(2,i) == 0 
            ID(i) = counter; %populating 'ID' array
            counter = counter + 1;
        end
    end
    Nun = max(ID);
    Nel = length(connections);
    
    if closed_boundaries(1)==1
        BC_D(end+1) = 1;
    elseif closed_boundaries(1)==2
        BC_S(end+1) = 1;
    else
        BC_N(end+1) = 1;
    end
    if closed_boundaries(2)==1
        BC_D(end+1) = Nnp;
    elseif closed_boundaries(2)==2
        BC_S(end+1) = Nnp;
    else
        BC_N(end+1) = Nnp;
    end
end