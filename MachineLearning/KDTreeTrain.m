function Tree = KDTreeTrain(Data,cols)

% Builds a KD tree to split Data.  The tree is grown by
% recursively splitting each node using the median value for that dimension
% until the leaf is down to one data point

% Data is an nxm matrix, where n is the number of points and m is the
% number of features.
% cols is a cell-vector of labels for each feature

% RETURNS t, a structure with three entries:
% t.p is a vector with the index of each node's parent node
% t.inds is the rows of X in each node (non-empty only for leaves)
% t.labels is a vector of labels showing the decision that was made to get
%     to that node
% t.Leaves: cell array 

% Create an empty decision tree, which has one node and everything in it
inds = {1:size(Data,1)}; % A cell per node containing indices of all data in that node
p = 0; % Vector contiaining the index of the parent node for each node
labels = {}; % A label for each node
Leaves = {}; %a struct to contain thresholds and next index in tree
depth = 1;
dims = size(Data,2);

% Create tree by splitting on the root
[inds, p, labels, Leaves] = KDTreeRecursive(Data, depth, dims, inds, p, labels, cols, 1, Leaves);


Tree.inds = inds;
Tree.p = p;
Tree.labels = labels;
Tree.Leaves = Leaves; %This is the actual tree used for online stuff