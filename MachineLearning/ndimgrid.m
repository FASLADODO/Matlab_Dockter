function Grid = ndimgrid(bounds,steps)
%create n-dimensional grid return coordinates of each point
%"bounds" is an 2xn matrix of the bounds for the grid
%rows = [min;max], columns = dimensions
%Use: bounds = DataBounds(X)
%"steps" is a 1xn vector of the numbers of steps in each dimension

%Returns an Mxn grid of points where M is prod(steps). Each row is a point
%in n-d space and comprises a complete grid
%Copyright Rodney Dockter 2017

%Get number of dimensions
[~,dims] = size(bounds);
[~,dss] = size(steps);
if(dss ~= dims)
    %if enough steps for each dimension is not given then duplicate first
    %dimension
   steps = steps(1).*ones(1,dims); 
end

%get evenly spaced points in each dimension
nsdim = max(steps); %max number of steps
dspace = zeros(nsdim,dims);
for dd = 1:dims
    dspace(:,dd) = linspace(bounds(1,dd),bounds(2,dd),nsdim)';
end

%compute total number of points
totalpoints = prod(steps);

%preallocate grid holder
Grid = zeros(totalpoints,dims);

%loop through all points in grid
for i = 1:totalpoints
    runningi=i-1;
    
    %all dimensions
    for j = 1:dims
       %idx = i/prod(steps(1:j-1)) mod steps(j) (Super efficient)
       Grid(i,j) = dspace( mod(runningi,steps(j) ) + 1 , j);
       runningi = floor(runningi / steps(j));
    end
    
end

%returns Grid
end