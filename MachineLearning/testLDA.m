% Note: discriminant coefficients are stored in W in the order of unique(Target)
%
% Example:
%

n1 = 20;
n2 = 24;
n3 = 17;

% Generate example data: 2 groups, of 10 and 15, respectively
X1 = randn(n1,2);
X2 =  randn(n2,2) + 2;
%X3 =  randn(n3,2) + 8;
X = [X1; X2];
Y = [ones(n1,1); ones(n2,1)*2];


%makin plots yo
figure
h1=scatter(X1(:,1),X1(:,2),'r+');
hold on
h2=scatter(X2(:,1),X2(:,2),'b*');
% hold on
% h3=scatter(X3(:,1),X3(:,2),'go');
hold off
xlabel('x1')
ylabel('x2')
legend([h1(1),h2(1)],{'class1','class2'})

% Calculate linear discriminant coefficients
[W,C] = LDASimple(X,Y);

%classify
[P,Class] = LDAonline(X,W,C);

%classification accuracy
classify = Class == Y;
accuracy = sum(classify)/length(classify)

%plot dat shiz
figure
scatter(X(:,1),X(:,2),20,P(:,1));
hold on
h = LDAplot(X,W,C);
hold off
colorbar
colormap cool
title('LDA boundary')