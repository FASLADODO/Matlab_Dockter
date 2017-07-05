function [match,Errors,T,ModelTransformed] = IterativeClosestPoint3D(DataModel,DataTarget)
%DataModel is 3D point set model
%DataTarget is 3D point set target
%http://ias.cs.tum.edu/_media/teaching/ss2012/techcogsys/as12perception-tutorial-solutions.pdf
    
%get data size
[NM,SM] = size(DataModel);
[NT,ST] = size(DataTarget);

% Will make Data1 and Dat2 be zero mean
muM = mean(DataModel);
muT = mean(DataTarget);
DataModelCenter = DataModel - repmat(muM,NM,1);
DataTargetCenter = DataTarget - repmat(muT,NT,1);

% compute matrix H to get rotation matrix
H = DataModelCenter' * DataTargetCenter;

%SVD to get rotation matrix
[U,S,V] = svd(H);
R = V*U';

%Use difference between centroids after rotation to get translation
t = muT' - R * muM';

%Now get transformation matrix by combining rotation and translation
T = [R, t; 0 0 0 1];

%now iterate
iter = 5;
for oo = 1:iter
    oo
    prevT = T;
    
    %Transform model into target data
    PadModel = [DataModel, ones(NM,1)]';
    transform = (T * PadModel)' ;
    DataModelTransform = transform(:,[1:3]);
    
    %KD tree
    kdOBJ = KDTreeSearcher(DataTarget);

    %Find nearest neighbors
    [match, mindist] = knnsearch(kdOBJ,DataModelTransform);
    
    %match rows
    Target2Model = DataTarget(match,:);
    
    %compute error
    Errors = NormRowWise(Target2Model - DataModelTransform);
    if(mean(Errors) < 0.01)
        iter
        break;
    end
    
    %figure out our new transformation matrix usign least squares
    PadModel = [DataModel, ones(NM,1)];
    for dd = 1:3
        outtemp = Target2Model(:,dd);
        paramrow = pinv(PadModel)*outtemp;
        T(dd,:) = paramrow';
    end
    %last row
    T(4,:) = [0,0,0,1];
    %scale rotation vectors
    for dd = 1:3
        T([1:3],dd) = T([1:3],dd) / sum(T([1:3],dd));
    end

end
ModelTransformed = (T*[DataModel,ones(NM,1)]')';

end