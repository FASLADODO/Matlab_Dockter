% Trying out Ahmidi's AFF-DCC method

r = 2;
t = -pi/2:0.1:pi/2;

x = r*cos(t)';
y = r*sin(t)';
z = t';

Data = [x,y,z];

figure
scatter3(Data(:,1),Data(:,2),Data(:,3),'ro')