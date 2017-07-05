function [scale,Grid,p] = integralRBF(Data,bw)
%integral estimate using gridded area and RBF

%get data size
[N,S] = size(Data);

%data size = number of steps in each dimension
nn = round(N^(1/S));

%data bounds in each dimensions
bounds = DataBounds(Data);

%make grid
Grid = ndimgrid(bounds,ones(1,S)*nn);

%get average step size of data
stepsize = diff(bounds)/nn;

%get probability at each point in grid
p = rbfpdist2(Grid,Data,bw);


%get the scale from the nd area times the height
scale = prod(stepsize)*sum(p);

if(false)
   figure
   scatter(Grid(:,1),p,'g.'); 
end

end