% Compute the ordinary Euclidean distance
D_act = 1

X = randn(100, 2);
Y = randn(100, 2) + D_act;
scatter(X(:,1),X(:,2),'rx')
hold on
scatter(Y(:,1),Y(:,2),'bx')
hold off
xlabel('x1')
ylabel('x2')

[IDX,DNN] = knnsearch(X,Y);

avg_D = mean(DNN)

%%

D = pdist2(X,Y,'euclidean'); % euclidean distance