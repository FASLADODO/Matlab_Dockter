function [inds, p, labels, Leaves] = KDTreeRecursive(Data,depth,dims,inds, p, labels, cols, node, Leaves)
% Recursively splits nodes based on median
% Data: data matrix with rows as samples
% depth: current dimensions
% dims: total dimensions
% inds: indices of all data in a given node
% p: parent node for each node
% labels: a label for each node
% cols: labels for each feature
% node: current node
% Leaves: struct to hold info for online
    depth
    %get data left on this branch
    current_Data = Data(inds{node},:);
    [NN,~] = size(current_Data);
    
    %check if we should leave
    if NN == 1
        Leaves{node}.Resultant = current_Data;
        return;
    end
    
    %figure out which dimension we are splitting on at this node
    axis = mod(depth-1,dims)+1;
    
    %temporary data
    dtemp = current_Data(:,axis);
    medval = median(dtemp);
    
    %find before and after data
    bid = dtemp < medval;
    aid = ~bid;
   
    inds = [inds; inds{node}(bid); inds{node}(aid)];
    inds{node} = [];
    
    %stash the current nodes parent and labels
    p = [p; node; node];
    labels = [labels; sprintf('%s < %2.2f', cols{axis}, medval); ...
    sprintf('%s >= %2.2f', cols{axis}, medval)];

    %put the thresholds and next index in Tree
    Leaves{node}.Threshold = medval;
    Leaves{node}.Feature = axis;
    Leaves{node}.FeatureLabel = cols{axis};
    Leaves{node}.Resultant = []; %leave empty if we don't classify on this node
    
    % Recurse on newly-create nodes
    n = numel(p)-2;
    Leaves{node}.ltindex = n+1; %store which index to go to
    [inds, p, labels, Leaves] = KDTreeRecursive(Data, depth+1, dims, inds, p, labels, cols, n+1, Leaves);
    Leaves{node}.gtindex = n+2; %store which index to go to
    [inds, p, labels, Leaves] = KDTreeRecursive(Data, depth+1, dims, inds, p, labels, cols, n+2, Leaves);
    
end