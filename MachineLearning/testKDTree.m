% test kd tree


Data = [2,3;5,4;9,6;4,7;8,1;7,2];

figure(10)
scatter(Data(:,1),Data(:,2))

%%
cols = {'d1', 'd2'};
t = KDTreeTrain(Data,cols)

%% Display the tree
figure
treeplot(t.p');
title('Decision tree ("**" is an inconsistent node)');
[xs,ys,h,s] = treelayout(t.p');

for i = 2:numel(t.p)
	% Get my coordinate
	my_x = xs(i);
	my_y = ys(i);

	% Get parent coordinate
	parent_x = xs(t.p(i));
	parent_y = ys(t.p(i));

	% Calculate weight coordinate (midpoint)
	mid_x = (my_x + parent_x)/2;
	mid_y = (my_y + parent_y)/2;

    % Edge label
	text(mid_x,mid_y,t.labels{i-1});
    
    % Leaf label
    if ~isempty(t.inds{i})
        val = Data(t.inds{i});
        if(isempty(t.Leaves{i}.Resultant))
            text(my_x, my_y, sprintf('y=%2.2f', t.Leaves{i}.Threshold));
        else
            text(my_x, my_y, sprintf('x=%2.2f, y=%2.2f', t.Leaves{i}.Resultant(1), t.Leaves{i}.Resultant(2)));
        end
    end
end