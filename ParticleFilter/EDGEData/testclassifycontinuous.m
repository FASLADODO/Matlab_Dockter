%% try classifying with pdfs

nn = 1000;

MU1 = [0 1];
SIGMA1 = [2 0; 0 .5];

MU2 = [-2 -1];
SIGMA2 = [1 0; 0 3];


X1 = mvnrnd(MU1,SIGMA1,nn);
X2 = mvnrnd(MU2,SIGMA2,nn);

cc = 1;
Data{cc}.State = X1;
Data{cc}.Class = ones(1,nn)*cc;
cc = 2;
Data{cc}.State = X2;
Data{cc}.Class = ones(1,nn)*cc;


figure
cc = 1;
scatter(Data{cc}.State(:,1),Data{cc}.State(:,2),10,'r.')
hold on
cc = 2;
scatter(Data{cc}.State(:,1),Data{cc}.State(:,2),10,'b.')
hold off

title('Training Distributions')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Group1','Group2','Group3')



%get gaussian mixture model
numgauss = 1;
cc = 1;
Models{cc}.GM = fitgmdist(Data{cc}.State,numgauss);
numgauss = 1;
cc = 2;
Models{cc}.GM = fitgmdist(Data{cc}.State,numgauss);


dlim = [];
for ii = 1:length(Data)
    dlim = [dlim; Data{ii}.State];
end
limz = [min(dlim(:,1)) max(dlim(:,1)) min(dlim(:,2)) max(dlim(:,2))];


NP = 50;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xplot = [XG1,XG2];

cc = 1;
Distribution{cc} = pdf(Models{cc}.GM,Xplot);
Models{cc}.Scale = 1/max(Distribution{cc});
Distribution{cc} = Distribution{cc} ./ max(Distribution{cc});


cc = 2;
Distribution{cc} = pdf(Models{cc}.GM,Xplot);
Models{cc}.Scale = 1/max(Distribution{cc});
Distribution{cc} = Distribution{cc} ./ max(Distribution{cc});


az = 20;
el = 30;
figure
cc = 1;
scatter(Data{cc}.State(:,1),Data{cc}.State(:,2),2,'r.');
hold on
m1 = Surface3Dalt(Xplot(:,1),Xplot(:,2),Distribution{cc},limz,NP);
set(m1,'facecolor','none')

hold on
cc = 2;
scatter(Data{cc}.State(:,1),Data{cc}.State(:,2),2,'g.');
hold on
m2 = Surface3Dalt(Xplot(:,1),Xplot(:,2),Distribution{cc},limz,NP);
set(m2,'facecolor','none')
hold off

colormap(cool);
colorbar;
title('Fit, Class 2')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);


%% Classifying

checkem = [];
numclasses = length(Data);

for ii = 1:numclasses
    for jj = 1:numclasses
        onlinePDF{ii}(:,jj) = pdf(Models{jj}.GM,Data{ii}.State) .* Models{jj}.Scale;
    end
    

    %HOW DO WE FIND OPTIMAL THRESHOLD (right now its hurestic)
    thresh = 10;
    for tt = 1:length(onlinePDF{ii})
        dda = onlinePDF{ii}(tt,:);
        [val,idx] = max( dda );
        temp = repmat(dda(:,idx),1,numclasses-1);
        others = dda;
        others(:,idx) = [];
        ratio = temp./others;
        check = ratio > thresh;
        Classify{ii}(tt,:) = [ii, check.*idx];
        Ratios{ii}(tt,:) = ratio;
        
        if(idx ~= ii && check)
            checkem = [checkem; ratio];
        end
    end
end

mean(checkem)


for ii = 1:numclasses
    temp = Classify{ii}(:,1) == Classify{ii}(:,2);
    Accuracy(ii) = sum(temp)/length(temp);
end
Accuracy

