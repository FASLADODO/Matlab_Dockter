function h = plotEllipse(xc,yc,rx,ry,n)
%plot ellipse

    theta = linspace(0,2*pi,n);
    x = rx*cos(theta) + xc;
    y = ry*sin(theta) + yc;

    h = plot(x,y,'k--','LineWidth',2);
end