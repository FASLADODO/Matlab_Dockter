%test Radial Basis Function

%From
%http://stats.stackexchange.com/questions/63881/use-gaussian-rbf-kernel-for-mapping-of-2d-data-to-3d

r1 = 2;
r2 = 4;
nn = 100;
nz = 0.3;

theta = linspace(0,2*pi,nn)';

d1 = [r1*sin(theta) + randn(nn,1)*nz, r1*cos(theta) + randn(nn,1)*nz];
d2 = [r2*sin(theta) + randn(nn,1)*nz, r2*cos(theta) + randn(nn,1)*nz];

figure
scatter(d1(:,1),d1(:,2),'r*');
hold on
scatter(d2(:,1),d2(:,2),'b+');
hold off
title('Original 2D data')

%now compute radial basis
gamma = 2;
z1 = RadialBasisFunction(d1,gamma);
z2 = RadialBasisFunction(d2,gamma);

%plot in 3D
figure
scatter3(d1(:,1),d1(:,2),z1,'r*');
hold on
scatter3(d2(:,1),d2(:,2),z2,'b+');
hold off
title('3D radial basis')


%% Try matlabs svm RBF
%http://www.mathworks.com/help/stats/fitcsvm.html

X = [d1;d2];
Y = [ones(nn,1)*(-1); ones(nn,1)*1]; 

SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');

[est,score] = predict(SVMModel,X);

plot(est)

Bounds = DataBounds(X);
ng = 50;
xl1 = linspace(Bounds(1,1),Bounds(2,1),ng);
xl2 = linspace(Bounds(1,2),Bounds(2,2),ng);
[XG1, XG2]  = meshgrid(xl1, xl2);


f = zeros(length(xl1), length(xl2));

gamma = SVMModel.KernelParameters.Scale

% Iter. all SVs
for i=1:length(SVMModel.SupportVectors)
    alpha_i = SVMModel.Alpha(i,:);
    sv_i    = SVMModel.SupportVectors(i,:);
    y_i     = SVMModel.SupportVectorLabels(i,:);
    for j=1:length(xl1)
        for k=1:length(xl2)
            x = [xl1(j),xl2(k)];
            f(j,k) = f(j,k) + alpha_i*y_i*constructKernel(x, sv_i, gamma, 'Gaussian');
            
        end
    end    
end


ZR = RadialBasisFunction(X,gamma);

figure
handle = Surface3D(X(:,1), X(:,2), ZR);


for i = 1:size(XG1, 2)
   this_X = [XG1(:, i), XG2(:, i)];
   [est,score] = predict(SVMModel,this_X);
   vals(:, i) = est;
end

%plot in 3D
figure
scatter(d1(:,1),d1(:,2),'r*');
hold on
scatter(d2(:,1),d2(:,2),'b+');
hold on
surf(XG1,XG2,f);
hold on
contour(XG1, XG2, vals, [0 1], 'Color', 'b');
hold off
title('3D radial basis')


figure
handle = plotSVMboundary(X,SVMModel);
% gscatter(XG1(:),XG2(:),vals(:));



