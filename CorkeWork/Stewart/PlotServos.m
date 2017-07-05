function PlotServos(theta,L,A1,A2)

if( length(theta) ~= 3 )
    disp('Too many angles')
end

J1 = [0,50,100;
         0,0,0];

dx = A1*cos(theta);
dy = -A1*sin(theta);

J2 = J1 + [dx;dy];

J3 = J1 + [ 0,0,0; L];

figure(2)
for ii = 1:3
   plot([J1(1,ii), J2(1,ii)] , [J1(2,ii), J2(2,ii)], 'b' ,'LineWidth',3 )
   hold on
   plot([J2(1,ii), J3(1,ii)] , [J2(2,ii), J3(2,ii)], 'r' ,'LineWidth',3 )
end

hold off
grid on

grid on
axis([-A1/2,max(max(J1))+(3*A1/2) ,-A1, A1+A2])
title('STU - Servo Plots')
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')


end