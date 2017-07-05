function [GridLimits] = SegmentLimitGrid(Grid)
% helper function to compute min and max bounds of a given grid element in
% n-dimensions
%Copyright Rodney Dockter 2017
[NN,SS] = size(Grid);

%get all possible permutations
steps = NN-1;
bounds = [ones(1,SS); ones(1,SS)*steps];
idxperm = round(ndimgrid(bounds,steps));

GridLimits = [];
for ii = 1:size(idxperm,1)
   %grab current permutation and the next in grid for max
   rowmin = idxperm(ii,:);
   rowmax = rowmin + 1;
   %grab that regions min and max
   regionmin = diag(Grid(rowmin,:));
   regionmax = diag(Grid(rowmax,:));
   GridLimits(ii,1,:) = regionmin;
   GridLimits(ii,2,:) = regionmax;
end

end