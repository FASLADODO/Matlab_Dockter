function [IV_Effort,Endpoints] = IntentVectorEffort(Data)
%compute the overall direction vector for the segment
%then compute the distance from that vector at each time step
%use this as measure of effort ie the farther from the intent the more
%effort

 %data size
[NN,SS] = size(Data);

%get starting and ending points
%get point1;
X1 = Data(1,:);
%get point2
X2 = Data(end,:);
Endpoints = [X1;X2];

%ultimate objective
IntentVector = UnitVec(X2-X1);

IV_Effort = zeros(NN,1);
%loop through all positions
for pp = 1:NN
    dtemp = Data(pp,:);
    %compute distance between that point and intent vector
    dist = 0; %TODO how to compute distance from point to vector in 3D
    %stash it
    IV_Effort(pp) = dist;
end

end