% http://www.mathworks.com/help/matlab/examples/graphs-and-matrices.html

[B,V] = bucky;
G = graph(B);

%plot the fig
figure
p = plot(G);
% p.XData = V(:,1);
% p.YData = V(:,2);

A = full(B);

figure
imshow(A)

%% Try with rand data

%make the data
rng default
X1 = randn(5,2);
X2 = randn(5,2) + repmat([2,2],5,1);

X = [X1; X2];

figure
scatter(X(:,1),X(:,2));



%construct adjacency matrix using k nearest neighbors 
[IDX,DST] = knnsearch(X,X,'k',3);

[NN,SS] = size(X);

A = zeros(NN,NN);
for ii = 1:NN-1
    %weights used to be 1, now are inverse of distances
    A(ii,IDX(ii,2)) = 1/DST(ii,2);
    A(IDX(ii,2),ii) = 1/DST(ii,2);
    A(ii,IDX(ii,3)) = 1/DST(ii,3);
    A(IDX(ii,3),ii) = 1/DST(ii,3);
end

%make it sparse
Asparse = sparse(A);

%plot the sparseness
figure
spy(Asparse)

%make graph from adjacency matrix
G = graph(Asparse);

%see the graph data
G.Edges

%plot the graph
figure
p = plot(G, 'EdgeLabel',G.Edges.Weight);
p.XData = X(:,1);
p.YData = X(:,2);






