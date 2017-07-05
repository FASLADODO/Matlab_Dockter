function [Progress] = IntentVectorProgress(Data)
%This gives a scalar for each data point which indicates the percent
%progress along the intent vector

%find the endpoints
E1 = Data(1,:);
E2 = Data(end,:);

%get the intent vector
Intent_Vec = E2 - E1;
Intent_Mag = norm(Intent_Vec);

%make this a unit vector
Intent_Unit = UnitVec(Intent_Vec);

%get data location all relative to first end point
Data_Relative = bsxfun(@minus,Data,E1);

%get dot product for each relative vector to the intent vector
Completion = Data_Relative*Intent_Unit';

%scale the completion by the overall distance
Progress = Completion/Intent_Mag;

end