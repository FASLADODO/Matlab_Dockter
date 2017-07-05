
%%
% hue 3d plots point matching

% read full image
picname = 'dummy2.jpg';
inputpic = imread(picname);
gray = rgb2gray(inputpic);
BWim = edge(gray,'sobel');
figure, imshow(BWim), title('BWim');
[row,col] = find(BWim);
impos = [row,col];
[rowsize,colsize,channels]=size(inputpic);

%read in model from edge image
picmodel = 'model3.jpg';
model = imread(picmodel);
graymodel = rgb2gray(model);
BWmod = im2bw(graymodel, 0.5);
figure, imshow(BWmod), title('BWmod');
[rowm,colm] = find(BWmod);
modelpos = [rowm,colm];

%%
% use icp algorithm
max_iter=400;
min_iter=40;
fitting=2; %[2,w]
thres=1e-5;
init_flag=1;
tes_flag=1;

[Rfit, Tfit] = icp(impos, modelpos, max_iter,min_iter,fitting,thres,init_flag,tes_flag);

Rfit
Tfit

T2 = Rfit'*Tfit % In original coordinates

%%
%draw circle of Tfit position of model image

for j = 1: size(rowm)
    newmodel(j,:) = (Rfit * modelpos(j,:)')';
end

plotrow = newmodel(:,1) + T2(1,1);
plotcol = newmodel(:,2) + T2(2,1);
figure, imshow(inputpic), title('onim');

hold on
plot(plotcol,plotrow,'s','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',2);
            
%%


























%