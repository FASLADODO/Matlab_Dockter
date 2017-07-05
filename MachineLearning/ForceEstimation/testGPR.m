%Test GPR for force estimation

[input,target,input_true,target_true] = createTrainingData();

[NN,SS] = size(input);

%plot the true vs the estimated
figure
scatter(1:NN,target,'bo')
hold on
scatter(1:NN,target_true,'r.')
hold off
legend('observed','true')


figure
scatter3(input(:,1),input(:,2),input(:,3),10,target)
xlabel('angle')
ylabel('velocty')
zlabel('acc')
colormap cool
h = colorbar;
ylabel(h, 'Torue Out');
title('sim data noisy')

figure
scatter3(input_true(:,1),input_true(:,2),input_true(:,3),10,target_true)
xlabel('angle')
ylabel('velocty')
zlabel('acc')
colormap cool
h = colorbar;
ylabel(h, 'Torue Out');
title('sim data true')


%% big subplots

figure
ax1 = subplot(5,1,1);
plot(input(:,1))
ylabel('angle')

ax2 = subplot(5,1,2); 
plot(input(:,2))
ylabel('vel')

ax3 = subplot(5,1,3); 
plot(input(:,3))
ylabel('acc')

ax4 = subplot(5,1,4); 
plot(input(:,4))
ylabel('T in')

ax5 = subplot(5,1,5); 
plot(target)
ylabel('T out')
suptitle('sim data noisey')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')


figure
ax1 = subplot(5,1,1);
plot(input_true(:,1))
ylabel('angle')

ax2 = subplot(5,1,2); 
plot(input_true(:,2))
ylabel('vel')

ax3 = subplot(5,1,3); 
plot(input_true(:,3))
ylabel('acc')

ax4 = subplot(5,1,4); 
plot(input_true(:,4))
ylabel('T in')

ax5 = subplot(5,1,5); 
plot(target_true)
ylabel('T out')
suptitle('sim data true')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')


%% Now lets mean variance scale the shit out of eveything

[X,xshift,xscale] = MeanVarianceScale(input_true);
[Y,yshift,yscale] = MeanVarianceScale(target_true);


figure
ax1 = subplot(5,1,1);
plot(X(:,1))
ylabel('angle')

ax2 = subplot(5,1,2); 
plot(X(:,2))
ylabel('vel')

ax3 = subplot(5,1,3); 
plot(X(:,3))
ylabel('acc')

ax4 = subplot(5,1,4); 
plot(X(:,4))
ylabel('T in')

ax5 = subplot(5,1,5); 
plot(Y)
ylabel('T out')
suptitle('scaled data')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

%% Train a GPR

coluse = [1,2,3,4]; %angle,vel,acc,torque
XTrain = X(:,coluse);
YTrain = Y;

sigma0 = std(YTrain);
sigmaF0 = mean(std(XTrain));
ld = 1500;
parmas = [ld,sigmaF0]

% gprMdl = fitrgp(X,Y,'KernelFunction','squaredexponential','KernelParameters',parmas,'Sigma',sigma0);
gprMdl = fitrgp(XTrain,YTrain,'KernelFunction','squaredexponential');

%% try it online

Xtest = X(:,coluse);
Ytest = predict(gprMdl,Xtest);

%get back into real units
% Yh = (Ytest ./ repmat(yscale,NN,1)) + repmat(yshift,NN,1);

figure
scatter(1:NN,YTrain,'bo')
hold on
scatter(1:NN,Ytest,'r.')
hold off
legend('observed','est')