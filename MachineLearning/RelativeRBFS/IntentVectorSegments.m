function [IV_Angle,Endpoints,IV_Angle_Diff] = IntentVectorSegments(Data)
%compute the overall direction vector for the segment
%then compute the angle of each individual motion to determine how close
%to the original motion we are at each time step

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

%steps in xyz
StepTemp = diff(Data);
StepVector = UnitVector(StepTemp); %make in unit

IV_Angle = zeros(NN-1,1);
%loop through all pos step vectors
for pp = 1:NN-1
    svec_temp = StepVector(pp,:);
    %compute angle between that vector and intent
    [angle,~,~] = AxisAngle3D(IntentVector,svec_temp);
    %stash it
    IV_Angle(pp) = angle;
end
%stash the intent vector angle
IV_Angle(end+1) = IV_Angle(end);

%stash the intent vector angle derivative
IV_Angle_Diff = diff(IV_Angle);
IV_Angle_Diff(end+1) = IV_Angle_Diff(end);

end