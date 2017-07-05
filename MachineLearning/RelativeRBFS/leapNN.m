
%% REQUIRES NN TOOLBOX (USE THE CITRIX MATLAB)

% Test leap motion data for seperability
% 
% D1 = load('2016-8-19_15-32-42_circles.csv');
% D2 = load('2016-8-19_15-34-47_ymotion.csv');

% D1 = load('2016-8-22_15-6-11_expert.csv');
% D2 = load('2016-8-22_15-6-57_novice.csv');

% D1 = load('2016-8-24_12-58-53_confident.csv');
% D2 = load('2016-8-24_13-5-35_timid.csv');

D1 = load('2016-8-30_12-53-20_confident.csv');
D2 = load('2016-8-30_12-53-52_timid.csv');


Data = [D1;D2];

Labels = [ones(length(D1),1)*1;ones(length(D2),1)*2];
% LabelKey = {'circles';'ymotion'};
LabelKey = {'expert';'novice'};


key.col.id = 1;
key.col.time = 2;
key.col.x = 3;
key.col.y = 4;
key.col.z = 5;
key.col.segment = 6;
key.col.dx = 7;
key.col.dy = 8;
key.col.dz = 9;
key.col.velmag = 10;
key.col.velalpha = 11;
key.col.velbeta = 12;
key.col.ddx = 13;
key.col.ddy = 14;
key.col.ddz = 15;
key.col.accmag = 16;
key.col.accalpha = 17;
key.col.accbeta = 18;

key.strings = {'frameid','time','x','y','z','segment','dx','dy','dz','velmag','velalpha','velbeta','ddx','ddy','ddz','accmag','accalpha','accbeta'};


% figure
% gscatter3(Data(:,key.col.x),Data(:,key.col.y),Data(:,key.col.z),Labels)
% title('initial data')


cslist = unique(Labels);

DataCopy = Data;

for cc = 1:length(cslist)
    %get current class data
    idx = find(Labels == cslist(cc) );
    Dtemp = DataCopy(idx,:);
    
    %start time from zero;
    Data(idx,key.col.time) = Data(idx,key.col.time) - Dtemp(1,key.col.time);
    
    %get time step
    dt = mean(diff(Dtemp(:,key.col.time)))
    
    %filter the motion just a lil bit
    xtemp = smooth(Dtemp(:,key.col.x));
    ytemp = smooth(Dtemp(:,key.col.y));
    ztemp = smooth(Dtemp(:,key.col.z));
    
    %compute derivatives
    dx = Calculate_velocity( xtemp, dt, 'holoborodko'); 
    dy = Calculate_velocity( ytemp, dt, 'holoborodko'); 
    dz = Calculate_velocity( ztemp, dt, 'holoborodko');
    
    ddx = Calculate_velocity( dx, dt, 'holoborodko'); 
    ddy = Calculate_velocity( dy, dt, 'holoborodko'); 
    ddz = Calculate_velocity( dz, dt, 'holoborodko');
    
    %magnitude and angle of velocity
    velmag = NormRowWise([dx,dy,dz]);
    [velalpha,velbeta] = PitchYaw3D([dx,dy,dz]);
    
    %magnitude and angle of acceleration
    accmag = NormRowWise([ddx,ddy,ddz]);
    [accalpha,accbeta] = PitchYaw3D([ddx,ddy,ddz]);
    
    %store that shiz back in the main matrix
    Data(idx,key.col.x) = xtemp;
    Data(idx,key.col.y) = ytemp;
    Data(idx,key.col.z) = ztemp;
    Data(idx,key.col.dx) = dx;
    Data(idx,key.col.dy) = dy;
    Data(idx,key.col.dz) = dz;
    Data(idx,key.col.velmag) = velmag;
    Data(idx,key.col.velalpha) = velalpha;
    Data(idx,key.col.velbeta) = velbeta;
    Data(idx,key.col.ddx) = ddx;
    Data(idx,key.col.ddy) = ddy;
    Data(idx,key.col.ddz) = ddz;
    Data(idx,key.col.accmag) = accmag;
    Data(idx,key.col.accalpha) = accalpha;
    Data(idx,key.col.accbeta) = accbeta;
    
    ClassData{cc} = Data(idx,:);
    
    totalSegments(cc) = max(Data(idx,key.col.segment)) - 1;
end

totalSegments

SegData = [];
%get segment data
for cc = 1:length(cslist)
    for ii = 1:totalSegments(cc)
        idx = find(ClassData{cc}(:,key.col.segment) == ii );
        SegData{cc}.segid{ii} = ClassData{cc}(idx,:);
    end
end

testcolumns = [key.col.dx,key.col.dy,key.col.dz,key.col.velalpha,key.col.velbeta,key.col.velmag, key.col.accmag,key.col.ddx,key.col.ddy,key.col.ddz,key.col.accalpha,key.col.accbeta,key.col.accmag];

%% Do a neural network using NN TOOLBOX 

X = Data(:,testcolumns)';
%get NN labels
targets = getNNLabels(Labels)';

%makes the nnet object
net = patternnet(20);
%view(net)

%train and view
net = train(net,X,targets);

%estimate classes
netout = net(X);
perf = perform(net,targets,netout);

%check accuracy
[~,classest] = max(netout);
corr = classest == Labels';
acc = mean(corr)

figure
plot(netout(1,:),'r')
hold on
plot(netout(2,:),'b')
hold off
legend('class 1','class 2')

%% now try classifying each segment
segclassify = [];

for cc = 1:length(cslist)
    for ii = 1:totalSegments(cc)
        tempd = SegData{cc}.segid{ii}(:,testcolumns)';
        netest = net(tempd);
        totalest = sum(netest,2);
        [~,classest] = max(totalest);
        segclassify = [segclassify; cc,classest];
    end
end




