function [center] = circlefit(y4,y1,y3,y2)
d=79.557;
n=25;
X=[n,n+d,n+d,n];
Y=[y1,y2,y3,y4];
xbar = (1/numel(X))*sum(X);
ybar = (1/numel(Y))*sum(Y);
u = X-xbar;
v = Y-ybar;
u2 = sum(u.*u);
uv = sum(u.*v);
u3 = sum(u.*u.*u);
uv2 = sum(u.*v.*v);
v2 = sum(v.*v);
v3 = sum(v.*v.*v);
vu2 = sum(v.*u.*u);
uc = ((1/2)*(v3+vu2)-((1/2)*(u3+uv2)*v2)/(uv))/(uv-((u2*v2)/uv));
vc = ((1/2)*(v3+vu2)-uc*uv)/(v2);
xc = vpa(uc+xbar);
yc = vpa(vc+ybar);
radius=vpa(sqrt(uc^2+vc^2+(u2+v2)/4));
center = [xc, yc, radius];




    