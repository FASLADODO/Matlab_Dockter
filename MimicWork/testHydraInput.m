d = load ('pick_and_place_2015_6_3_16_27_2.txt');

[steps,inputs] = size(d)

%%

%time
t = d(:,1);
dt = diff(t);

if( inputs >= 7)
    %positions
    X1 = d(:,2);
    Y1 = d(:,3);
    Z1 = d(:,4);
    X2 = d(:,5);
    Y2 = d(:,6);
    Z2 = d(:,7);
end

if(inputs >= 15)
    %orientations (quaternion)
    q1_0 = d(:,8);
    q1_1 = d(:,9);
    q1_2 = d(:,10);
    q1_3 = d(:,11);

    q2_0 = d(:,12);
    q2_1 = d(:,13);
    q2_2 = d(:,14);
    q2_3 = d(:,15);
    
    EA1 = SpinCalc('QtoEA123',[q1_0,q1_1,q1_2,q1_3],0.05,1) .* (pi/180);
    EA2 = SpinCalc('QtoEA123',[q2_0,q2_1,q2_2,q2_3],0.05,1) .* (pi/180);
end



figure(1)

h1 = quiver3(X1,Y1,Z1,EA1(:,1),EA1(:,2),EA1(:,3),2,'r');
hold on
h2 = quiver3(X2,Y2,Z2,EA2(:,1),EA2(:,2),EA2(:,3),2,'k');

hold off
grid on
title('segmented tool motion')
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
legend([h1, h2],'Tool1 grasping (right)','Tool2 grasping (left)')