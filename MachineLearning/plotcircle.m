function h = plotcircle(center,r,Color)
%Plot circle on 2D Plot
%center = [X,Y]
%r = Radius scalar
%Color = optional color argument
    if(nargin < 3)
       Color = 'black'; 
    end
    x = center(1);
    y = center(2);
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'Color',Color);
end