%%test sort

X = randn(100,2);
scatter(X(:,1),X(:,2))

scaleX = X - repmat(mean(X),length(X),1);

mean(scaleX)

%%
XD = NormRowWise(X);

tic 
[XDS,idxs] = sort(XD);
toc 

XR = X(idxs,:);

scatter3(XR(:,1),XR(:,2),XDS);