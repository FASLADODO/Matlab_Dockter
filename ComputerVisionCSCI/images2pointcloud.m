% Create a 3D point cloud from a series of images

filenames = {'SquareTemplate1.png'; 'SquareTemplate2.png'; 'SquareTemplate3.png'};

NF = length(filenames);

for ff = 1:NF
    filenames{ff}
    
    RGB = imread(filenames{ff});
    input{ff}.mat = rgb2gray(RGB);
    figure, imshow(input{ff}.mat)
end


%% now get true data points in x,y,z pointclouds

AllPoints = [];
zscale = 1.5;
zgranulairty = 10; % 10 z points per layer

%%loop through each image
for ff = 1:NF
    [NR,NC] = size(input{ff}.mat); %rows and columns
    
    for yy = 1:NR
        for xx = 1:NC
            val = input{ff}.mat(yy,xx);
            if(val  > 125)
                zheight = (ff-1)*zscale;
                zpoints = linspace(zheight,zheight+zscale,zgranulairty);
                for zz = zpoints;
                    tempd = [xx,yy,zz,ff];
                    AllPoints = [AllPoints; tempd];
                end
            end
        end
    end
end

fsize=14;
labstr = {'Layer 1','Layer 2','Layer 3'}

figure
gscatter3(AllPoints(:,1),AllPoints(:,2),AllPoints(:,3),labstr(AllPoints(:,4)))
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)
xlabel('x','FontSize',fsize);
ylabel('y','FontSize',fsize);
zlabel('z','FontSize',fsize);

%% add noise to point cloud


shift = [0.5,0.5,0.5];
noise = 0.1;

[np,sp] = size(AllPoints);

AllPointsNoise = AllPoints(:,[1:3]) + randn(np,3).*noise + repmat(shift,np,1);

figure
gscatter3(AllPointsNoise(:,1),AllPointsNoise(:,2),AllPointsNoise(:,3),labstr(AllPoints(:,4)))
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)
xlabel('x','FontSize',fsize);
ylabel('y','FontSize',fsize);
zlabel('z','FontSize',fsize);

%% try ICP between the two

[Correspondence,Errors,T,transformdata] = IterativeClosestPoint3D(AllPoints,AllPointsNoise);
T

figure
scatter3(transformdata(:,1),transformdata(:,2),transformdata(:,3),'b+')
hold on
scatter3(AllPointsNoise(:,1),AllPointsNoise(:,2),AllPointsNoise(:,3),'ro')
hold off



