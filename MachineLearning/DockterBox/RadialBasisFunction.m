function Z = RadialBasisFunction(data,gamma)
%Calculate Gaussian Radial Basis Function
%data: Rows are Samples, columns are dimensions
%Gamma = 2

%From
%http://stats.stackexchange.com/questions/63881/use-gaussian-rbf-kernel-for-mapping-of-2d-data-to-3d

%pdist finds distance from each row of data to all other rows of data
%each row is considered a coordinate in N-D space
D = squareform( pdist(data, 'euclidean') );
D = exp(-(D .^ 2) ./ ( 2*gamma^2));
Z = sum(D)' / length(data); %scale by length yo

% D = NormRowWise(data);
% Z = exp(-(D .^ 2) ./ ( 2*gamma^2));

end

