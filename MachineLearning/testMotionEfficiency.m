%testing motion efficiency function

% SIMPLE LINE
nn = 100;
x = [ [1:nn]', [1:nn]' ];

M = [3;4];
B = 2;

yo = x*M + B;
z = yo + randn(nn,1)*10;

Data = [x,z];
[D,X1,X2,params] = MotionEfficieny(Data);

figure
scatter3(x(:,1),x(:,2),z)
hold on
scatter3(X1(:,1),X1(:,2),X1(:,3),'go')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),'go')
hold off
title('sample data')

figure
scatter3(x(:,1),x(:,2),z,20,D)
hold on
endp = [X1;X2];
plot3(endp(:,1),endp(:,2),endp(:,3),'g');
hold off
title('motion efficiency')
colorbar 
colormap cool


%%  weird sinusoid

nn = pi*30;
x = [ [0:nn]', [0:nn]' ];

M = [2;2];
B = 2;

%linear function
yo = x*M + B;
%add in a siusoid
z = yo + 20*sin(0.5*x(:,1));

Data = [x,z];
[D,X1,X2,params] = MotionEfficieny(Data);

figure
scatter3(x(:,1),x(:,2),z)
hold on
scatter3(X1(:,1),X1(:,2),X1(:,3),'go')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),'go')
hold off
title('sample data')

figure
scatter3(x(:,1),x(:,2),z,20,D)
hold on
endp = [X1;X2];
plot3(endp(:,1),endp(:,2),endp(:,3),'g');
hold off
title('motion efficiency')
colorbar 
colormap cool

%%  big arc

nn = 100;
x = [ [0:nn]', [0:nn]' ];

M = [2;2];
B = 0;
center = [200,200];
radius = 400;

%linear function
z = x*M + B;
%add in a circle to x and y
x(:,1) = radius*cos( (x(:,1)/nn)*1.57);
x(:,2) = radius*sin( (x(:,2)/nn)*1.57);


Data = [x,z];
[D,X1,X2,params] = MotionEfficieny(Data);

figure
scatter3(x(:,1),x(:,2),z)
hold on
scatter3(X1(:,1),X1(:,2),X1(:,3),'go')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),'go')
hold off
title('sample data')

figure
scatter3(x(:,1),x(:,2),z,20,D)
hold on
endp = [X1;X2];
plot3(endp(:,1),endp(:,2),endp(:,3),'g');
hold off
title('motion efficiency')
colorbar 
colormap cool

