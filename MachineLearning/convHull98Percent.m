function [DT,hull] = convHull98Percent(Data,Sigma)
%Data: ND matrix containing training data
%Sigma percent of data to get inside hull eg: 'Sigma = 0.98'
%DT is the deluanay triangulation
%hull: are the vertices of the hull

%Get the index length corresponding to 98%
limit = round(Sigma*length(Data));

%center the data around a mean
centered = Data - repmat(mean(Data),length(Data),1);

%get distances to center
distz = NormRowWise(centered);

%sort all points according to distance
[~,idx] = sort(distz);

%Get back data in sorted order
sdata = Data(idx,:);

%get limited data (98% of all data)
ldata = sdata(1:limit,:);

%comute the deluanay triangulation
DT = delaunayTriangulation(ldata); 

%get the convex hull from DT
[hull,~] = convexHull(DT);

end