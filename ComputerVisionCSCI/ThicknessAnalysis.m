% analyse layer via transparency

input = imread('IMG_2942.JPG');

figure, imshow(input)

%% make it grayscale

gray = rgb2gray(input);

figure, imshow(gray)

%% identify regions

layer1idx{1} = [200:1500];
layer1idx{2} = [100:630]; %becuz row column
layer2idx{1} = [200:1500];
layer2idx{2} = [630:1060];
layer3idx{1} = [200:1500];
layer3idx{2} = [1060:1450];

xcol = 2;
ycol = 1;
layer1im = gray(layer1idx{xcol}, layer1idx{ycol});
layer2im = gray(layer2idx{xcol}, layer2idx{ycol});
layer3im = gray(layer3idx{xcol}, layer3idx{ycol});

figure, imshow(layer1im)
figure, imshow(layer2im)
figure, imshow(layer3im)

%% measure the transparency in each region


layer1val = mean(mean(double(layer1im)))
layer1std = std(std(double(layer1im)))
layer2val = mean(mean(double(layer2im)))
layer2std = std(std(double(layer2im)))
layer3val = mean(mean(double(layer3im)))
layer3std = std(std(double(layer3im)))


%% make some dope ass box and whiskers

layer1all = double(layer1im(:));
layer2all = double(layer2im(:));
layer3all = double(layer3im(:));


fsize = 14;
figure
C = [layer1all;layer2all;layer3all];
grp = [ones(length(layer1all),1); ones(length(layer2all),1)*2; ones(length(layer3all),1)*3];
boxplot(C,grp,'Notch','on','Labels',{'Layer 1','Layer 2','Layer 3'})
ylabel('Transparency','FontSize',fsize)












