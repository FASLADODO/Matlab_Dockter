clear all
bounds = [45, 24, 1; 89, 70, 100];
steps = [10 , 10, 10];

tic
Grid = ndimgrid(bounds,steps);
toc

scatter3(Grid(:,1),Grid(:,2),Grid(:,3),'r.')