%Just SVM again (quadtratic programming)

%create some fake data
u1 = [2,3];
u2 = [3,5];
c1 = [0.5,0.5];
c2 = [0.5,0.5];
nn = 100;

X0 = [c1(1)*randn(nn,1) + u1(1), c1(2)*randn(nn,1) + u1(2)];
X1 = [c2(1)*randn(nn,1) + u2(1), c2(2)*randn(nn,1) + u2(2)];

%Get working data
X = [X0;X1];
Y = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)

%TRy the lagrangian dual
% THIS is thegeneral method
% alpha = zeros(length(Y),1);
% 
% sumtemp = 0;
% for ii = 1:NN
%     for jj = 1:NN
%         sumtemp = sumtemp + alpha(ii)*alpha(jj)*Y(ii)*Y(jj)*X(ii,:)*X(jj,:)' ;
%     end
% end
% 
% W = sum(alpha) - (1/2)*sumtemp;
% %constrain alpha >= 0
% %constrain sum(alpha.*Y) == 0;
%x_i with non-zero alpha_i are called support vectors (SV) 
%solve with quadratic programming


[W,B,SV] = SVM_train(X,Y)

xp = linspace(min(X(:,1)), max(X(:,1)), 100);
yp = - (W(1)*xp + B)/W(2); %THIS IS THE SHIT

%ezplot way
f = @(x,y) B + [x, y]*W';

figure
scatter(X0(:,1),X0(:,2),'r.');
hold on
scatter(X1(:,1),X1(:,2),'b.');
hold on;
h = ezplot(f, [min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))]);
% plot(xp, yp, '-g'); 
hold off
title('classify')

[Accuracy,Classify] = SVM_Online(X,Y,W,B);
Accuracy