function [LABELS] = NormalizedGraphCuts(Data,nn)
%perform spectral clustering using normalized graph cuts
%X= n dimensional data matrix with columns as dimensions
%nn = number of groupings

%see \MachineLearning\graphcuts.m for example usage
%http://www.sci.utah.edu/~gerig/CS7960-S2010/handouts/Normalized%20Graph%20cuts.pdf
% http://stackoverflow.com/questions/15480635/spectral-clustering

%now this is our adjacency matrix for the optimal cut
dist2 = pdist2(Data,Data);
W = exp(-dist2);
sumw = sum(W,2);
D = diag(sumw);

% (D-W)x=lambda Dx 
%solve this as A*V = B*V*D using [V,D]=eig(A,B)
L=D-W;
L=D^(-0.5)*L*D^(-0.5);
[eVEC,eVAL]=eig(L);

V1=eVEC(:,1:nn) % choose the k eigenvectors

%So graph cuts actually just makes k means search easier
[~,LABELS] = kMeansIterative(V1,nn,100);

end