function integral_result = IntegralND(func,bounds,nsteps)
%integral estimate using gridded reimann sum estimate
%func is an anonymous function evaluated column wise for each dimension
%bounds is a column wise min and max in each dimension 
%nsteps = number of steps in each dimension of grid

% eg for 2d:
% func = @(x) exp(- sqrt( x(:,1)^2 + x(:,2).^2 ) ); % exp(- sqrt(x^2+y^2))
% bounds = [-1,-1;1,1]; % -1<x<1, -1<y<1
% nsteps = [10,10]; 10 steps in x, 10 steps in y

%get data size
[~,S] = size(bounds);

%set number of steps same in each dimension if not specified
if(length(nsteps) ~= S)
   nsteps = ones(1,S)*nsteps(1); 
end

%make grid
Grid = ndimgrid(bounds,nsteps);

%get average step size of data
stepsize = diff(bounds)./nsteps;

%evaluate the function at each point in the grid
evalf_grid = func(Grid);

%get the integral from the nd area times the height
integral_result = prod(stepsize)*sum(evalf_grid);

%for test plotting
if(false)
    figure
    scatter(Grid(:,1),Grid(:,2),'ro')
    hold on
    Surface3D(Grid(:,1),Grid(:,2),evalf_grid);
    hold off
    
end

end