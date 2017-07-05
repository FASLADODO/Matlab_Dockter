function [pt1,pt2] = circle_intersect(x1,y1,r1,x2,y2,r2)

A = [x1 y1]; %# center of the first circle
B = [x2 y2]; %# center of the second circle
a = r1; %# radius of the SECOND circle
b = r2; %# radius of the FIRST circle
c = norm(A-B); %# distance between circles

cosAlpha = (b^2+c^2-a^2)/(2*b*c);

u_AB = (B - A)/c; %# unit vector from first to second center
pu_AB = [u_AB(2), -u_AB(1)]; %# perpendicular vector to unit vector

%# use the cosine of alpha to calculate the length of the
%# vector along and perpendicular to AB that leads to the
%# intersection point
intersect_1 = A + u_AB * (b*cosAlpha) + pu_AB * (b*sqrt(1-cosAlpha^2));
intersect_2 = A + u_AB * (b*cosAlpha) - pu_AB * (b*sqrt(1-cosAlpha^2));

pt1 = intersect_1;
pt2 = intersect_2;

end
