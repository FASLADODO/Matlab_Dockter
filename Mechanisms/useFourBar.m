% Enter lengths

a = 2;
b = 3;
c = 5;
d = 3;
theta2 = 130;


hf = figure('color','white');
axis([-5 15 -5 15])

for j = 1:60
    theta2 = theta2-1;
    [x1,y1,x2,y2,x3,y3,x4,y4] = Four_Bar_Rod(a,b,c,d,theta2);
    line([x1,x2],[y1,y2],'LineWidth',4,'Color',[0 0 0])
    hold on
    line([x2,x3],[y2,y3],'LineWidth',4,'Color',[0 0 1])
    hold on
    line([x3,x4],[y3,y4],'LineWidth',4,'Color',[0 1 0])
    hold on
    line([x1,x4],[y1,y4],'LineWidth',4,'Color',[1 0 0])
    hold off

    refreshdata(hf,'caller')
    drawnow
    pause(1);
    clf
end