function h = convHull3DPlot(DT,hull)
%plots a 3D convex hull
h = trisurf(hull,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),...
       'FaceAlpha',0.1,'FaceColor','g');
end