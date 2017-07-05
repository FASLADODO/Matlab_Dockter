
nn = 1000;
D = abs(rand(nn,2));

metric0 = min(D,[],2);
metric1 = mean(D,2)*3 - abs(D(:,1) - D(:,2));
metric2 = mean(D,2)*4 - abs(D(:,1) - D(:,2));
metric3 = mean(D,2)*8 - abs(D(:,1) - D(:,2));
metric4 = mean(D,2);

figure
scatter(D(:,1),D(:,2),'r.')
hold on
Surface3D(D(:,1),D(:,2),metric0,'contour')
title('metric0')

figure
scatter(D(:,1),D(:,2),'r.')
hold on
Surface3D(D(:,1),D(:,2),metric1,'contour')
title('metric1')

figure
scatter(D(:,1),D(:,2),'r.')
hold on
Surface3D(D(:,1),D(:,2),metric2,'contour')
title('metric2')

figure
scatter(D(:,1),D(:,2),'r.')
hold on
Surface3D(D(:,1),D(:,2),metric3,'contour')
title('metric3')

figure
scatter(D(:,1),D(:,2),'r.')
hold on
Surface3D(D(:,1),D(:,2),metric4,'contour')
title('metric4')


