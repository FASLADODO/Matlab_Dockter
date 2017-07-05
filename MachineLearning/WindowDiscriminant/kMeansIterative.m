function [means,grpIndex] = kMeansIterative(X,km,maxiter)
%kmeans using iterative voronoi groupings
%Inputs:
%X= data matrix with each column as a state, each row is observations
%km= the number of means to search for ie km=3
%maxiter= optional arguments for maximum iterations
%Outputs:
%means= mean coordinates which will be [km x #states]
%grpIndex= corresponding mean index for each observation in X


%From https://sites.google.com/site/dataclusteringalgorithms/k-means-clustering-algorithm

    if(nargin == 2)
       maxiter = 1000; 
    end

    [NN,SS] = size(X);
    
    %initialize means
    rng('shuffle')
    seeds = randi(NN, km, 1);
    means = X(seeds,:);

    %initialize groupings
    grpIndex = ones(NN,1);
    prevgroups = grpIndex;
    dist2means = zeros(NN,km);
    
    %checks for convergence
    doneso = 0;
    
    %go through all iterations
    for tt = 1:maxiter

        %compute voronoi groups
        for mm = 1:km
            scale = X - repmat(means(mm,:),NN,1);
            dist = sqrt(sum(abs(scale).^2,2));

            dist2means(:,mm) = dist;
        end

        %get inidices for minimums (determines center groupings)
        [~,grpIndex] = min(dist2means,[],2);


        %use data groupings to compute new means
        for mm = 1:km
            groupdata = X(grpIndex == mm, :);

            means(mm,:) = mean(groupdata,1);
        end

        %check if we stopped changing
        if(isequal(grpIndex,prevgroups) )
           doneso = 1;
           fprintf('convergence in %d steps \n',tt)
           break; 
        end

        prevgroups = grpIndex;
    end


    if(doneso == 0)
       disp('failed to converge!') 
    end

end