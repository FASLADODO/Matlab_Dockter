function classification = BoundingBoxCheck(Data,limits)
    %Checks if set of points in Data are inside the cube defined by limits
    %returns 1 for inside and outside the box

    %Data: data matrix with columns as states
    %limits: 2xn matrix of lower and upper bounds of cube in n dimensions
    % (limits comes from PlotBoundingCube)

    [NN,SS] = size(Data);
    [NL,SL] = size(limits);
    if(SS ~= SL)
       error('columns in Data must match columns in limit') 
    end

    %check if our data is below the min or above the max in any dimensions
    within = [Data <= repmat(limits(1,:),NN,1), Data >= repmat(limits(2,:),NN,1) ];

    %sum up check in each rows
    alldims = sum(within,2);
    alldims(alldims > 0) = 1;

    %final labels
    classification = alldims;
    classification = ~classification; 

end