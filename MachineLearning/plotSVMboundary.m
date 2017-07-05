function handle = plotSVMboundary(Data,Model)
% Given a model from fitcsvm() and a Data matrix this will plot the boundary colors

steps = 100;

[NN,SS] = size(Data);

% Get bounds and grid
bounds = DataBounds(Data);
Grid = ndimgrid(bounds,repmat(steps,1,SS));

[NG,SG] = size(Grid);

for i = 1:NG
   this_X = Grid(i, :);
   [est,score] = predict(Model,this_X);
   vals(i,:) = est;
end

handle = [];

gid = unique(vals);
colors = 'rgb';
markers = 'osd';

if(SS == 2)
    handle = gscatter(Grid(:,1), Grid(:,2), vals);
elseif(SS == 3)
    for idx = 1 : length(gid)
        data = Grid(vals == gid(idx),:);
        handle = plot3(data(:,1), data(:,2), data(:,3), [colors(idx) markers(idx)]);
        hold on;
    end
end


end