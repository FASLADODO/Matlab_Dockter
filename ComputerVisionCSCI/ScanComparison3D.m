% load up and compare 3D scan

data = load('hand_pyramid_moving.xyzrgb'); %still scan
%data = load('hand_pyramid_still.xyzrgb'); %still scan

figure
scatter3(data(:,1),data(:,2),data(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
title('original scan')
axis equal

%% do the transformation stuff

PosData = data(:,[1:3]);
[NN,SS] = size(PosData);


%Build rotation matrices
theta = deg2rad(24); %30 for still
alpha = deg2rad(-3); % 0
RX = [1,0,0,0;
    0,cos(theta),-sin(theta),0;
    0,sin(theta),cos(theta),0;
    0,0,0,1];
RY = [cos(alpha),0,sin(alpha),0;
    0,1,0,0;
    -sin(alpha),0,cos(alpha),0;
    0,0,0,1];
T = RX*RY;


%transform the data
PadData = [PosData, ones(NN,1) ];
TransformData = (T*PadData')';
TransformData = TransformData(:,[1:3]);

%get min
minheight = min(TransformData(:,3));
TransformData(:,3) = TransformData(:,3) -minheight;

% peak
[maxheight,id] = max(TransformData(:,3));
center = TransformData(id,:);

figure
scatter3(TransformData(:,1),TransformData(:,2),TransformData(:,3),5,TransformData(:,3))
hold on
scatter3(center(:,1),center(:,2),center(:,3),'r+')
title('transformed scan')
xlabel('x')
ylabel('y')
zlabel('z')



%% shift data

ShiftData = TransformData;
ShiftData(:,[1:2]) = TransformData(:,[1:2]) - repmat(center(:,[1:2]),NN,1);

fsize = 16;
figure
scatter3(ShiftData(:,1),ShiftData(:,2),ShiftData(:,3),5,ShiftData(:,3))
%title('centered scan')
xlabel('x (mm)','FontSize',fsize)
ylabel('y (mm)','FontSize',fsize)
zlabel('z (mm)','FontSize',fsize)


%% floor it to get rid of johns weird flour ring

boxlim = 32;

xlimplus = 20;
xlimminus = -25;
ylimplus = 20;
ylimminus = -20;

%find all our data in the hump and out
% idin = find( sqrt(ShiftData(:,1).^2 + ShiftData(:,2).^2 ) <= boxlim);
% idout = find( sqrt(ShiftData(:,1).^2 + ShiftData(:,2).^2 ) > boxlim);
idin = find( ShiftData(:,1) > xlimminus & ShiftData(:,1) < xlimplus & ShiftData(:,2) > ylimminus & ShiftData(:,2) < ylimplus );
idout = find( ShiftData(:,1) < xlimminus | ShiftData(:,1) > xlimplus | ShiftData(:,2) < ylimminus | ShiftData(:,2) > ylimplus );

%grab it
datain = ShiftData(idin,:);
dataout = ShiftData(idout,:);

%limit it
ShiftData(idout,3) = 0;
ShiftData(idin,3) = ShiftData(idin,3) - min(ShiftData(idin,3));

fsize = 16;
figure
scatter3(ShiftData(:,1),ShiftData(:,2),ShiftData(:,3),5,ShiftData(:,3))
%title('centered scan')
xlabel('x (mm)','FontSize',fsize)
ylabel('y (mm)','FontSize',fsize)
zlabel('z (mm)','FontSize',fsize)

%% interpolation

scaler = 23.62; % 600 dpi
step = 1 / scaler; % steps per mm

[Xq,Yq] = meshgrid(-50:step:50);

F = scatteredInterpolant(ShiftData(:,1),ShiftData(:,2),ShiftData(:,3));
Vq = F(Xq,Yq);

InterpData = [Xq(:),Yq(:),Vq(:)];

figure
mesh(Xq,Yq,Vq);
% title('interpolated scan')
xlabel('x (mm)','FontSize',fsize)
ylabel('y (mm)','FontSize',fsize)
zlabel('z (mm)','FontSize',fsize)
% axis equal

[maxheight,id] = max(InterpData(:,3));

%% also pull our brim up

% InterpData(InterpData(:,3) < 0,3) = 0;

figure
scatter3(InterpData(:,1),InterpData(:,2),InterpData(:,3),5,InterpData(:,3))
% title('interpolated scan')
xlabel('x (mm)','FontSize',fsize)
ylabel('y (mm)','FontSize',fsize)
zlabel('z (mm)','FontSize',fsize)

%% now load up our images (MODEL)


% Create a 3D point cloud from a series of images

filenames = {'pyramid_layer1.png'; 'pyramid_layer3.png'; 'pyramid_layer5.png'; 'pyramid_layer7.png'; 'pyramid_layer9.png'};

NF = length(filenames);

for ff = 1:NF
    filenames{ff}
    
    RGB = imread(filenames{ff});
    input{ff}.mat = rgb2gray(RGB);
    figure, imshow(input{ff}.mat)
end

%% find min and max in each image

%%loop through each image
for ff = 1:NF
    [NR,NC] = size(input{ff}.mat); %rows and columns
    
    minx = NC;
    maxx = 0;
    miny = NR;
    maxy = 0;
    for yy = 1:NR
        for xx = 1:NC
            val = input{ff}.mat(yy,xx);
            if(val < 125)
                if(xx < minx)
                    minx = xx;
                end
                if(xx > maxx)
                    maxx = xx;
                end
                if(yy < miny)
                    miny = yy;
                end
                if(yy > maxy)
                    maxy = yy;
                end
            end
        end
    end
    %outer bounds
    OB{ff}.minx = minx;
    OB{ff}.maxx = maxx;
    OB{ff}.miny = miny;
    OB{ff}.maxy = maxy;
end


%% now get true data points in x,y,z pointclouds

totallayers = 15;
ModelPoints = [];
zscale = maxheight/ totallayers;

layercopy = 3;
layeridx = 0;

%%loop through each image
for ff = 1:NF
    [NR,NC] = size(input{ff}.mat); %rows and columns
    
    %3 layers per image
    for cc = 1:layercopy
        for yy = 1:NR
            for xx = 1:NC
                val = input{ff}.mat(yy,xx);
                if(val < 125 && (OB{ff}.minx == xx || OB{ff}.maxx == xx || OB{ff}.miny == yy || OB{ff}.maxy == yy) )
                    zheight = (layeridx)*zscale + zscale;
                    tempd = [xx,yy,zheight];
                    ModelPoints = [ModelPoints; tempd];
                elseif (val < 125 && layeridx == 14)
                    zheight = (layeridx)*zscale  + zscale;;
                    tempd = [xx,yy,zheight];
                    ModelPoints = [ModelPoints; tempd];
                end
            end
        end
        layeridx = layeridx + 1;
    end
end

fsize=14;

figure
scatter3(ModelPoints(:,1),ModelPoints(:,2),ModelPoints(:,3),5,ModelPoints(:,3))
xlabel('x','FontSize',fsize);
ylabel('y','FontSize',fsize);
zlabel('z','FontSize',fsize);


%% tack on some zeros for interp

[XZeros, YZeros] = meshgrid(0:2:100);
datzeros = [XZeros(:),YZeros(:)];
ZZeros = zeros(length(datzeros),1);
datzeros = [datzeros, ZZeros];

%remove middle
rmid = find(datzeros(:,1) >= OB{1}.minx & datzeros(:,1) <= OB{1}.maxx & datzeros(:,2) >= OB{1}.miny & datzeros(:,2) <= OB{1}.maxy);
datzeros(rmid,:) = [];

ModelPoints = [datzeros;ModelPoints];


figure
scatter3(ModelPoints(:,1),ModelPoints(:,2),ModelPoints(:,3))


%% center scale

[NM,SM] = size(ModelPoints);
centermodel = [50,50];

ShiftModel = ModelPoints;
ShiftModel(:,[1:2]) = ModelPoints(:,[1:2]) - repmat(centermodel,NM,1);

figure
scatter3(ShiftModel(:,1),ShiftModel(:,2),ShiftModel(:,3),5,ShiftModel(:,3))
title('centered model')
xlabel('x')
ylabel('y')
zlabel('z')



%% interpolation

scaler = 23.62; % 600 dpi
step = 1 / scaler; % steps per mm

[Xq,Yq] = meshgrid(-50:step:50);

F = scatteredInterpolant(ShiftModel(:,1),ShiftModel(:,2),ShiftModel(:,3));
Vq = F(Xq,Yq);

InterpModel = [Xq(:),Yq(:),Vq(:)];

fsize = 16
figure
mesh(Xq,Yq,Vq);
%title('interpolated model')
xlabel('x (mm)','FontSize',fsize)
ylabel('y (mm)','FontSize',fsize)
zlabel('z (mm)','FontSize',fsize)


%% compare the two
figure
h1 = scatter3(InterpData(:,1),InterpData(:,2),InterpData(:,3),'r.')
hold on
h2 = scatter3(InterpModel(:,1),InterpModel(:,2),InterpModel(:,3),'b.')
hold off
legend([h1(1),h2(1)],'Data','Model')
xlabel('x')
ylabel('y')
zlabel('z')



%% fine tune

[NMH,SMH] = size(InterpModel);
[NDH,SDH] = size(InterpData);

shiftfine = [2,2,0];

ModelSub = InterpModel;
DataSub = InterpData + repmat(shiftfine,NDH,1);

figure
h1 = scatter3(DataSub(:,1),DataSub(:,2),DataSub(:,3),'r.')
hold on
h2 = scatter3(ModelSub(:,1),ModelSub(:,2),ModelSub(:,3),'b.')
hold off
legend([h1(1),h2(1)],'Data','Model')
xlabel('x')
ylabel('y')
zlabel('z')


%% run volumetric differences

[NM,SM] = size(ModelSub);
ColumnArea = step*step

%demoninators for stuff
TotalDeposit = 0; 
TotalModel = 0; 

% For error heat map
ErrorMap = zeros(NM,1);

TPSum = 0; %Deposited where we we wanted
FPSum = 0; %deposited more than we wanted
TNSum = 0; %didnt deposit where we wanted
for ii = 1:NM
    tempData = DataSub(ii,3);
    tempModel = ModelSub(ii,3);
    
    diff= 0;
    if(DataSub(ii,1) < 20 && DataSub(ii,1) > -20 && DataSub(ii,2) < 30 && DataSub(ii,2) > -20)
        %Account for TP, TN, FP
        if(tempData > tempModel)
            diff = tempData - tempModel;
            TPSum = TPSum + tempModel; %inked where we wanted
            FPSum = FPSum + diff; %extra where we didnt want
            
        elseif(tempModel > tempData)
            diff = tempModel - tempData;
            TPSum = TPSum + tempData; %inked sorta where we wanted
            TNSum = TNSum + diff; %didnt ink all that we needed
        end
        %increment totals
        TotalDeposit = TotalDeposit + tempData;
        TotalModel = TotalModel + tempModel;
        
        %Error placement
        ErrorMap(ii,:) = diff;
    end
end

TotalDeposit = TotalDeposit*ColumnArea;
TotalModel = TotalModel*ColumnArea;
TPSum = TPSum*ColumnArea;
TNSum = TNSum*ColumnArea;
FPSum = FPSum*ColumnArea;

TP = TPSum / TotalModel
TN = TNSum / TotalModel
FP = FPSum / TotalDeposit

%% Error heat map

fsize = 16;
figure
Surface3D(DataSub(:,1),DataSub(:,2),ErrorMap,'contour');
xlabel('x (mm)','FontSize',fsize)
ylabel('y (mm)','FontSize',fsize)
colormap spring
hc = colorbar
ylabel(hc, 'Error (mm)','FontSize',fsize)

%% find match

%KD tree
kdOBJ = KDTreeSearcher(DataSub);

%Find nearest neighbors
[match, mindist] = knnsearch(kdOBJ,ModelSub);

%match rows
Data2Model = DataSub(match,:);

Errors = NormRowWise(Data2Model - ModelSub);

mean(Errors)