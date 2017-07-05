%With linear functions

clear all


%%
nn = 1000;
noiz = 0.5;

rangen = [1 30];

params1 = [1.5 25]; %m, b
params2 = [2.7 10]; %m, b


xm = linspace(rangen(1),rangen(2),nn)';

xm1 = xm + randn(nn,1)*noiz;
xm2 = xm + randn(nn,1)*noiz;

yx1 = xm1.*params1(1) + params1(2) + randn(nn,1)*noiz;
yx2 = xm2.*params2(1) + params2(2) + randn(nn,1)*noiz;


X1 = [xm1, yx1];
X2 = [xm2, yx2];

Data = [X1;X2];
Y = [ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
classes = unique(Y);

figure
scatter(X1(:,1),X1(:,2),'r.')
hold on
scatter(X2(:,1),X2(:,2),'b.')
%hold on
%scatter(X3(:,1),X3(:,2),'g.')
hold off
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);


std(X1)
std(X2)

(range(Data)./length(Data)) .* std(Data)

mean(sqrt(std(Data)))
    
%% plot the resulting seperations

[h,probgrid] = windowDiscriminantAlt(Data, Y);


%% Train the models for each
pthresh = 0.4;
kmeanlimit = 6;

if(~exist('kmg'))
    kmg = [];
end
[Models, kmg] = windowDiscriminantTrain(Data, Y, pthresh, kmg, kmeanlimit, 'ploton');

%% Log Likelihood ratios

disp('classification accuracy with Rods Discriminant:')
LLR_X1 = windowDiscriminantOnline(X1, Models, classes, pthresh);
LLR_X1(end);
class_okay1 = LLR_X1 >= 0;
Accuracy1 = sum(class_okay1)/length(class_okay1)

LLR_X2 = windowDiscriminantOnline(X2, Models, classes, pthresh);
LLR_X2(end);
class_okay2 = LLR_X2 <= 0;
Accuracy2 = sum(class_okay2)/length(class_okay2)

%% least squares accuracy

disp('classification accuracy with gaussians:')
AccLC = LeastSqauresClassify( [Data(:,1),ones(length(Data),1)], Data(:,2), Y);

figure
plot(AccLC(:,1),'r');
figure
plot(AccLC(:,2),'b');


figure
subplot(2,1,1);
plot([1:length(LLR_X1)],LLR_X1,'r-*' )
ylabel('Kernel LLR','FontSize', 12)
title('Discriminant Window vs. Least Squares Class 1','FontSize', 12)

subplot(2,1,2);
plot([1:length(AccLC(:,1))],AccLC(:,1),'b--+' )
ylabel('Least Squares','FontSize', 12)



figure
subplot(2,1,1);
plot([1:length(LLR_X2)],LLR_X2,'r-*' )
ylabel('Kernel LLR','FontSize', 12)
title('Discriminant Window vs. Least Squares Class 2','FontSize', 12)

subplot(2,1,2);
plot([1:length(AccLC(:,2))],AccLC(:,2),'b--+' )
ylabel('Least Squares','FontSize', 12)

