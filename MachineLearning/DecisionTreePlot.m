function [] = DecisionTreePlot(Tree,Labels)
% Display the tree
treeplot(Tree.p');
title('Decision tree');
[xs,ys,h,s] = treelayout(Tree.p');

for i = 2:numel(Tree.p)
	% Get my coordinate
	my_x = xs(i);
	my_y = ys(i);

	% Get parent coordinate
	parent_x = xs(Tree.p(i));
	parent_y = ys(Tree.p(i));

	% Calculate weight coordinate (midpoint)
	mid_x = (my_x + parent_x)/2;
	mid_y = (my_y + parent_y)/2;

    % Edge label
	text(mid_x,mid_y,Tree.labels{i-1});
    
    % Leaf label
    if ~isempty(Tree.inds{i})
        val = Labels(Tree.inds{i});
        if numel(unique(val))==1
            text(my_x, my_y, sprintf('y=%2.2f\nn=%d', val(1), numel(val)));
        else
            %inconsistent data
            text(my_x, my_y, sprintf('**y=%2.2f\nn=%d', mode(val), numel(val)));
        end
    end
end


end