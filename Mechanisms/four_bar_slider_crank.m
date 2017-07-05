%Slider crank inversion 1


L2 = 10;
L3 = 20;


theta2 = 120;
theta3 = 0;

x1 = 0;
y1 = 0;

y3 = 0;

hf = figure('color','white');

step = -1;


while 1
    theta2 = theta2 + step;
    theta3 = asind((L2/L3)*sind(theta2));
    L4 = L2*cosd(theta2) + L3*cosd(theta3);
    
    x2 = L2*cosd(theta2);
    y2 = L2*sind(theta2);
    x3 = L4;
    
    line([x1,x2],[y1,y2],'LineWidth',4,'Color',[0 0 0])
    hold on
    line([x2,x3],[y2,y3],'LineWidth',4,'Color',[0 0 1])
    hold on
    line([x3-3,x3+3],[y3,y3],'LineWidth',6,'Color',[0 1 0])
    hold off
    axis([-10 30 -10 30])
    
    refreshdata(hf,'caller')
    drawnow
    pause(0.001);
    clf
    
    if(theta2 == 15)
        step = 1;
    end
    
    if(theta2 == 121)
       step  = -1; 
    end
end