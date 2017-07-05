function [inds, p, labels, Leaves] = split_node(X, Y, inds, p, labels, cols, node, Leaves)
%https://alliance.seas.upenn.edu/~cis520/wiki/index.php?n=Lectures.DecisionTrees
% Recursively splits nodes based on information gain
% X: data matrix with rows as samples
% Y: labels for each sample
% inds: indices of all data in a given node
% p: parent node for each node
% labels: a label for each node
% cols: labels for each feature
% node: current node
% Leaves: struct to hold info for online

% Check if the current leaf is consistent (if data is perfectly classified)
if numel(unique(Y(inds{node}))) == 1
    current_Y = Y(inds{node});
    Leaves{node}.Resultant = mode(current_Y);
    return;
end

% Check if all inputs have the same features
% We do this by seeing if there are multiple unique rows of X
if size(unique(X(inds{node},:),'rows'),1) == 1
    current_Y = Y(inds{node});
    Leaves{node}.Resultant = mode(current_Y);
    return;
end

% Otherwise, we need to split the current node on some feature

%initialize some variables to find best feature
best_ig = -inf; %best information gain
best_feature = 0; %best feature to split on
best_val = 0; % best value to split the best feature on
stepratio = 0.5;

%Get the data left on this branch for this node (some data has already been
%partitioned from previous nodes)
current_X = X(inds{node},:);
current_Y = Y(inds{node});
[NN,SS] = size(current_X);

% Loop over each feature to decide which feature will be used for this
for i = 1:SS
    %current feature
    feat = current_X(:,i);
    
    % Deterimine the values to split on
    vals = unique(feat);
    if numel(vals) < 2
        continue
    end
    %weak learner uses steps between each data point for thresholds
    splits = stepratio*(vals(1:end-1) + vals(2:end));
    
    % Try splitting at each of the thresholds
    bin_mat = double(repmat(feat, [1 numel(splits)]) < repmat(splits', [NN 1]));
    
    % Compute the information gains for each threshold
    H = EntropySample(current_Y);
    H_cond = zeros(1, size(bin_mat,2));
    for j = 1:size(bin_mat,2)
        %H_cond(j) = ConditionalEntropy(curr_Y, bin_mat(:,j));
        H_cond(j) = cond_ent(current_Y, bin_mat(:,j));
    end
    IG = H - H_cond;
    
    % Find the best threshold
    [val ind] = max(IG);
    if val > best_ig
        best_ig = val;
        best_feature = i;
        best_val = splits(ind);
    end
end

% Split the current node into two nodes using threshold
% get indices for each branch
datt = current_X(:,best_feature);
feat = datt < best_val;
notfeat = ~feat;

%stash this into indices data for next leaves
inds = [inds; inds{node}(feat); inds{node}(notfeat)];
inds{node} = [];

%stash the current nodes parent and labels
p = [p; node; node];
labels = [labels; sprintf('%s < %2.2f', cols{best_feature}, best_val); ...
    sprintf('%s >= %2.2f', cols{best_feature}, best_val)];

%put the thresholds and next index in Tree
Leaves{node}.Threshold = best_val;
Leaves{node}.Feature = best_feature;
Leaves{node}.FeatureLabel = cols{best_feature};
Leaves{node}.Resultant = []; %leave empty if we don't classify on this node

% Recurse on newly-create nodes
n = numel(p)-2;
Leaves{node}.ltindex = n+1; %store which index to go to
[inds, p, labels, Leaves] = split_node(X, Y, inds, p, labels, cols, n+1, Leaves);
Leaves{node}.gtindex = n+2; %store which index to go to
[inds, p, labels, Leaves] = split_node(X, Y, inds, p, labels, cols, n+2, Leaves);

