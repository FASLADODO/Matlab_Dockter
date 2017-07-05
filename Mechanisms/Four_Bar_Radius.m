%Four bar using interseciton of circles method

%Specify this stuff...
L2 = 3;
L3 = 5;
L4 = 7;
L1 = 6;
theta1 = 0;
theta2 = 150;
elbow = 0; %1 for up, 0 for down


x1 = 0;
y1 = 0;
x3 = 0;
y3 = 0;
x4 = L1*cosd(theta1);
y4 = L1*sind(theta1);
x2 = L2*cosd(theta2);
y2 = L2*sind(theta2);

[pt1,pt2] = circle_intersect(x2,y2,L3,x4,y4,L4);

hf = figure('color','white');

while 1
    
    theta2 = theta2 - 1;
    x2 = L2*cosd(theta2);
    y2 = L2*sind(theta2);
    [pt1,pt2] = circle_intersect(x2,y2,L3,x4,y4,L4);

    if(elbow == 1)
       x3 = pt1(1);
       y3 = pt1(2);
    else
       x3 = pt2(1);
       y3 = pt2(2);
    end
    
    theta3 = atand((y3 - y2)/(x3-x2));
    theta4 = atand((y3 - y4)/(x3-x4));
    
    line([x1,x2],[y1,y2],'LineWidth',4,'Color',[0 0 0])
    hold on
    line([x2,x3],[y2,y3],'LineWidth',4,'Color',[0 0 1])
    hold on
    line([x4,x3],[y4,y3],'LineWidth',6,'Color',[0 1 0])
    hold off
    axis([x1 - L2, x4 + L3, y4 - L2, y4 + L3])

    refreshdata(hf,'caller')
    drawnow
    pause(0.001);
    clf

end