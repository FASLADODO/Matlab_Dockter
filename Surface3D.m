function handle = Surface3D(x,y,z,option,limits)
%Surface3D Summary: plot 3D surface from x,y,z vectors

nn = 30;
%linear spacing
if(nargin == 5)
    xlin = linspace(limits(1,1),limits(1,2),nn);
    ylin = linspace(limits(2,1),limits(2,2),nn);
else
    xlin = linspace(min(x),max(x),nn);
    ylin = linspace(min(y),max(y),nn);
end

%mesh grid
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(x,y,z,'natural', 'none');
Z = f(X,Y);
Z(isnan(Z)) = 0 ;

%check args
if(nargin < 4)
   option = 'mesh'; 
end

%plot
if(strcmp(option, 'mesh') )
    handle = mesh(X,Y,Z);
elseif(strcmp(option, 'contour') )
    handle = contourf(X,Y,Z);
elseif(strcmp(option, 'surf') )
    handle = surf(X,Y,Z);
else
    handle = mesh(X,Y,Z);
end

if(~strcmp(option, 'contour') )
    set(handle,'facecolor','none')
end
colormap(cool);
colorbar;

end

