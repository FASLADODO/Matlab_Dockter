% create some rand data


nn = 100;
ss = 3;

X = [];
Y = [];

%create some random data
s1 = 1;
mu1 = [4,5,4];
s2 = 1;
mu2 = [1,1,1];
x1 = randn(nn,ss).*s1 + repmat(mu1,nn,1);
x2 = randn(nn,ss).*s2 + repmat(mu2,nn,1);

%combining
X = [x1;x2];
Y = [ones(nn,1)*-1; ones(nn,1)*1 ];

classes = unique(Y);
Ynum = MapValues(Y,classes,[1,2]);

%plot the original classes
figure
gscatter3(X(:,1),X(:,2),X(:,3),Ynum);
title('initial data')

%% built in function

% T = clusterdata(X,'maxclust',3);
T = clusterdata(X,'linkage','ward','savememory','on','maxclust',2);

figure
scatter3(X(:,1),X(:,2),X(:,3),10,T)

%% Try pdist2 version

%http://www.sci.utah.edu/~gerig/CS7960-S2010/handouts/Normalized%20Graph%20cuts.pdf
% http://stackoverflow.com/questions/15480635/spectral-clustering

nn = 10;
%make the data
X1 = randn(nn,2);
X2 = randn(nn,2) + repmat([3,3],nn,1);
X3 = randn(nn,2) + repmat([6,6],nn,1);

X = [X1; X2; X3];
Y = [ones(nn,1)*1;ones(nn,1)*2;ones(nn,1)*3];

figure
gscatter(X(:,1),X(:,2),Y);
title('true class')

%DO GRAPH CUT
tic
%now this is our adjacency matrix for the optimal cut
dist2 = pdist2(X,X);
W = exp(-dist2);
sumw = sum(W,2);
D = diag(sumw);

% (D-W)x=lambda Dx 
%solve this as A*V = B*V*D using [V,D]=eig(A,B)
L=D-W;
L=D^(-0.5)*L*D^(-0.5);
[eVEC,eVAL]=eig(L);

k = 3; %number of desired groups
V1=eVEC(:,1:k) % choose the k eigenvectors

%So graph cuts actually just makes k means search easier
[~,LABELS] = kMeansIterative(V1,k,100);
toc

LABYA = NormalizedGraphCuts(X,3);

%plot the graph
figure
gscatter(X(:,1),X(:,2),LABYA);
title('Est class')


% MAKE GRAPHS FOR FUN

%construct adjacency matrix using k nearest neighbors 
[IDX,DST] = knnsearch(X,X,'k',3);

[NN,SS] = size(X);

A = zeros(NN,NN);
for ii = 1:NN-1
    %weights used to be 1, now are inverse of distances
    A(ii,IDX(ii,2)) = exp(-DST(ii,2));
    A(IDX(ii,2),ii) = exp(-DST(ii,2));
    A(ii,IDX(ii,3)) = exp(-DST(ii,3));
    A(IDX(ii,3),ii) = exp(-DST(ii,3));
end


%make it sparse
Asparse = sparse(A);

%plot the sparseness
figure
spy(Asparse)
title('sparseness')

%make graph from adjacency matrix
G = graph(Asparse);

%see the graph data
G.Edges

%plot the graph
figure
p = plot(G, 'EdgeLabel',G.Edges.Weight);
p.XData = X(:,1);
p.YData = X(:,2);
title('graph adjacency')


%% Now look at inter class pdist info


nn = 10;
%make the data
X1 = randn(nn,2);
X2 = randn(nn,2) + repmat([3,3],nn,1);

X = [X1; X2];
Y = [ones(nn,1)*1;ones(nn,1)*2];

figure
gscatter(X(:,1),X(:,2),Y);
title('true class')

% get groups
classes = unique(Y);
D = [];
for i = 1:length(classes)
    D{i} = X(Y == classes(i),:);
end

% get inter and intra
sep = [];
for i = 1:length(classes)
    to = pdist2(D{i},D{i});
    sep{i}.Inter = exp(-to);
    for j = 1:length(classes)
        if i ~= j
            td = pdist2(D{i},D{j});
            sep{i}.Intra{j} = exp(-td);
        end
    end
end

%TODO USE THIS FOR SOMETHING
dist2 = pdist2(X,X);
W = exp(-dist2);
sumw = sum(W,2);
D = diag(sumw);

% (D-W)x=lambda Dx 
%solve this as A*V = B*V*D using [V,D]=eig(A,B)
L=D-W;
L=D^(-0.5)*L*D^(-0.5);
[eVEC,eVAL]=eig(L);

k = 2; %number of desired groups
V1=eVEC(:,1:k) % choose the k eigenvectors

%So graph cuts actually just makes k means search easier
[~,LABELS] = kMeansIterative(V1,k,100);

%plot the graph
figure
gscatter(X(:,1),X(:,2),LABELS);
title('est class')



%% bio informatics

% http://www.mathworks.com/help/bioinfo/examples/working-with-graph-theory-functions.html

load oscillatorgraph

whos g names

spy(g)

gObj = biograph(g,names);

go = view(gObj);

[S,C] = conncomp(gObj);