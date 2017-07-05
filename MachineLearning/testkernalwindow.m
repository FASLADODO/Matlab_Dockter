nn = 100;

X = rand(nn,2)*20;



%% Test kernel window
rind = nn / 2;

%Args
CPoint =  X(rind,:);
windowSize = 7;

windowData = DataInWindow(X,CPoint,windowSize)


figure
scatter(X(:,1),X(:,2),'b.')
hold on
scatter(windowData(:,1),windowData(:,2),'r.')
hold on
scatter(CPoint(:,1),CPoint(:,2),'g*')
hold on
plotcircle(CPoint,windowSize);
hold off
axis equal
title('Window Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('Training Data','Kernel Data','Current Data Point','Window Bounds')
set(hl,'FontSize',12);


%% Now test window probs

for ii = 1:nn
   datatemp = X(ii,:); 
   [Pt] = probabilityWindow(X,datatemp);
   All_Prob(ii,1) = Pt;
end

figure
scatter(X(:,1),X(:,2),'b.')
hold on
handle = Surface3D(X(:,1),X(:,2),All_Prob)
hold off


