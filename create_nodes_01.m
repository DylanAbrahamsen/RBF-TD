function xt = create_nodes_01 

% Output parameter
%   xt  All nodes in the x,t-plane, sorted for increasing t-values

% addpath('G:\Research\RBF\Node_distributions\Matlab_04_trui\')
% addpath('G:\Research\RBF\RBF-FD_Time_stepping\Matlab_01\')

% Create the scattered node set in the x,t-plane
ninit  = 1e6;
dotmax = 1e7;
xt = node_drop ([0 1 0 1], ninit, dotmax, @radius_variable_02);

% Create boundary nodes
corners = [0 0;1 0;1 1;0 1];        % The four corners
bdyr = discretize_bdy (corners,@radius_variable_02);
bdy = bdyr(:,1:2);
xt2  = repel(xt,bdy,corners,@radius_variable_02); % Repel by all boundaries
xt = xt2;   % Remove last column with radius information

[ni,~] = size(xt);      % Get total number of interior nodes
[nb,~] = size(bdy);     % Get total number of boundary nodes

xt = [bdy;xt];          % Make xy contain both boundary and interior nodes
xt = sortrows(xt,2);    % Sort the rows from lowest y (i.e. t) to highest