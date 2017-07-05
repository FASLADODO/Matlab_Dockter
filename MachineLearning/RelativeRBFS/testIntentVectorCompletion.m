
%make some data

nn = 100;

x = [[1:nn]', [1:nn]'];

M = [0.2;0.2];

z0 = x*M;
z = z0 + randn(nn,1)*2;
Data = [x,z];

figure
scatter3(Data(:,1),Data(:,2),Data(:,3),'ro')
axis([0 100 0 100 0 100])

Progress = IntentVectorProgress(Data);

figure
scatter3(Data(:,1),Data(:,2),Data(:,3),10,Progress)
axis([0 100 0 100 0 100])
colorbar 
colormap cool

%% rotation matrix way (slow but sorta works)

%compute unit vector
primary = Data(end,:) - Data(1,:);
primvec = UnitVec(primary);
pitch = atan2(primvec(3),primvec(1));
yaw = atan2(primvec(2),primvec(1));

rad2deg(pitch)
rad2deg(yaw)

RY = [cos(yaw),-sin(yaw),0;sin(yaw),cos(yaw),0;0,0,1]
RP = [cos(pitch),0,sin(pitch);0,1,0;-sin(pitch),0,cos(pitch)]

mapData1 = RP*Data';
mapData = RY*mapData1;
mapData = mapData';
mapData1 = mapData1';

figure
scatter3(mapData1(:,1),mapData1(:,2),mapData1(:,3),'ro')

figure
scatter3(mapData(:,1),mapData(:,2),mapData(:,3),'ro')


%% weird way
%compute unit vector
primary = Data(end,:) - Data(1,:);
primvec = UnitVec(primary)
xaxis = [1,0,0];

%compute parameters
v = cross(primvec,xaxis)
s = norm(v)
c = dot(primvec,xaxis)
VX = [0,-v(3),v(2);;v(3),0,-v(1);-v(2),v(1),0]
R = eye(3) + VX + VX.*((1-c)/(s^2))
% R = [1,0,0;0,c,-s;0,s,c];

mapData = R*Data';
mapData = mapData';

figure
scatter3(mapData(:,1),mapData(:,2),mapData(:,3),'ro')

