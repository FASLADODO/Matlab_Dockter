function P = circle3D(r,C,O,U,ploton)
%C = [x;y;z] (center of circle)
%N = [alpha,beta,gamma] (rool pitch yaw)
%U is a unit vector from center to a point on the circle

%http://mathforum.org/library/drmath/view/63755.html

%Example Use:
% C = [5;5;5];
% O = [0.1,0.1,0.1];
% r = 3;
% P = circle3D(r,C,O,1)

%plot steps
steps = 0.01;

%create normal vector
N = [cos(O(3))*sin(O(2));
    sin(O(3))*sin(O(2));
	cos(O(2)) ];

%array for loop
t = 0:steps:2*pi;

%cross product
V = cross(N,U);

%create vector P for plotting
for ii = 1:length(t)
    P(:,ii) = r*cos(t(ii))*U + r*sin(t(ii))*V + C;
end

%for plotting
if(ploton)
    plot3(P(1,:),P(2,:),P(3,:),'r','LineWidth',3);
end

end