% Creates mesh properties. Specific for this project's geometry.
% - input 'h' is width and height of a single element in meters. preferably should be a simple fraction.
% - input 'width' is width of membrane in meters.
% - input 'height' is height of membrane in meters.
% - output 'nodes' is array of coordinates (x,y) of each node, and a closed/open code.
% - output 'connections' is array of corresponding nodes for each element. Identical to IEN array.
% - output 'ID' is array of open variables.
% - output 'FE_nodes' is array of nodes on F-E side (for boundary conditions).
% - output 'ED_nodes' is array of nodes on E-D side (for boundary conditions). (not used)
% - output 'DC_nodes' is array of nodes on D-C side (for boundary conditions). (not used)

function [nodes,connections,ID] = meshGenerator(h, width, height)
    global Nnp nodes_per_col F_node D_node F_position D_position;
    Nnp = (width / h + 1)*(height / h + 1);
    Nel = (width / h )*(height / h ) * 2;
    nodes = zeros(3, Nnp); % initialize
    connections = zeros(3, Nel); % initialize
    nodes_per_row = width/h + 1;
    nodes_per_col = height/h + 1;
    ID = zeros(1,Nnp);

    i = 1;
    for col = 1 : 1 : nodes_per_row
        for row = 1 : 1 : nodes_per_col
            nodes(1,i) = (col-1) * h; %populating 'nodes' array
            nodes(2,i) = (row-1) * h;
            nodes(3,i) = row == 1 || row == nodes_per_col || col == 1 || col == nodes_per_row; % code for open/closed node
            i = i + 1;
        end
    end
    
    i = 1;
    for node=1:1:Nnp  %populating 'connections' array
        if nodes(1,node) < width && nodes(2,node) < height
            % 'top' shape elements
            connections(1,i) = node;
            connections(2,i) = node + nodes_per_col + 1;
            connections(3,i) = node + 1;
            i = i + 1;
            % 'bottom' shape elements
            connections(1,i) = node;
            connections(2,i) = node + nodes_per_col;
            connections(3,i) = node + nodes_per_col + 1;
            i = i + 1;
        end
    end
    
    counter = 1;
    for i = 1:Nnp
        if nodes(3,i) == 0 
            ID(i) = counter; %populating 'ID' array
            counter = counter + 1;
        end
    end
end