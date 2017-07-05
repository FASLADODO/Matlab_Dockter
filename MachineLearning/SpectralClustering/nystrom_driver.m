function labels = nystrom_driver(data,sigma,nsamp,k,n,t)
% This function calculates the nystomized approximation to the
% Normalized cut partitioning algorithm.
%
% data is the matrix with the data, where each row of the matrix
% contains a data point
% 
% sigma is the standard deviation for the gaussian 
%
% nsamp is the number samples that should whie performing the
% approximation
% k is the number of desired groups
% n is the maximum number of iterations you would like the kmeans
% algorithm to run on top of the eigenvectors
%
% t is the threshold at which you would like the kmeans algorithm
% to terminate
%
% this function assumes the availability of a function
% "your_own_distance"
%
% which will take two matrices A and B, with row vectors
% corresponding to datapoints and calculate the all pairs distance
% between them.
%
% an example is included in this directory, its called "pairdist"
%
% labels is an array of cluster labels for each point.


[N,Ncol]=size(data); % the dimensions of the dataset
                   % each row is a new data point
                   
ind_all=randperm(N);
ind_samp=ind_all(1:nsamp);
ind_rest=ind_all(nsamp+1:end);
%samp_pos=[x(ind_samp)' y(ind_samp)'];



%d2A=your_own_distance(FB(ind_samp,:),FB(ind_samp,:)); % use your own distance function to get the pairwise distances
%d2B=your_own_distance(FB(ind_samp,:),FB(ind_rest,:));

%This is just pdist2(), why you write your own?
d2A=pairdist(data(ind_samp,:),data(ind_samp,:),'L2'); % use your own distance function to get the pairwise distances
d2B=pairdist(data(ind_samp,:),data(ind_rest,:),'L2');

%Weights nearby points higher, far away points lower
A=exp(-d2A/(2*sigma^2)); % choose your sigma 
B=exp(-d2B/(2*sigma^2));
%imagesc(A);
% the joys of nystrom

[V,ss]=nystrom_ncut(A,B);

% depermute the nsytrom result
V(ind_all,:)=V;
V1=V(:,1:k); % choose the k eigenvectors
Vnormalized=V1./repmat(sqrt(sum(V1.^2,2)),1,k);
[centers,labels,error]=km(Vnormalized',k,n,t);

