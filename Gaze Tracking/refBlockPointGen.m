%% Generates a set of points linearly interpolated between Ref row entries 
% at a rate of pointsPerUnitDist (e.g., points per cm or inches) 
% this means that if there is zero distance between two Ref points, no 
% points are added in that region; 
% Also, this means points should be spaced by distances of 
% 1/pointsPerUnitDist.  However, if consecutive points in Ref are spaced 
% less than that distance, they are stuffed in anyway, i.e. all points in 
% Ref shoud appear in data with some extra points thrown in at a rate of 
% pointsPerUnitDist.
%   Ref = [m x n] m consecutive n-dimensional points
%   pointsPerUnitDist = [double], how many points to add per cm or inch
%   data = [q x n] sequence of q n-dimensional points where q>=m;
%
%   data = refBlockPointGen(Ref, pointsPerUnitDist)
%
function data = refBlockPointGen(Ref, pointsPerUnitDist)

% Find lengths between points; (length to previous point)
deltaD = [ 0; sum(diff(Ref).^2,2).^.5 ];

% preallocate
data = zeros(floor(sum(deltaD)*pointsPerUnitDist), size(Ref, 2));

% set first point.
%data(1,:) = Ref(1,:);
countSoFar =0;

% Iterate through intermediate points
for i=1:size(Ref, 1)
    % number of points to add in this segment
    % if dcimal value larger than 0, round up
    segN = ceil(pointsPerUnitDist*deltaD(i)); 
    
    %segP=[];
    % for all points in the segment, find linear interp. for each coor pair 
    %for j = 1:segN
    
    % only process nonzero length segments (ignore redundant points)
    if(segN > 0)
        % for every column of Ref,k  ...
        for k = 1:size(Ref,2)
            p = linspace(Ref(i-1,k), Ref(i,k), segN+1); % add extra point since firs point will be discarded
            data(countSoFar+[1:(segN+1)], k) = p(1:(end)); %discard the first point
        end
        countSoFar = countSoFar + segN;
    end
    %data(i,:) =  
end


% set last point.
data(end, :) = Ref(end,:);