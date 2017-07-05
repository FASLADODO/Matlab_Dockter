function idx = getNewLandmark(Y)
    %Y all observation measurements
    bestdist = inf;
    thetabounds = 45 * (pi/180);
    idx = 1;
    for ii = 1:size(Y,2)
        %get angle and distance to current landmark
        theta2mark = Y(2,ii);
        thetaerror = abs(theta2mark);
        dist2mark = Y(1,ii);
        %check if this best meets our needs
        if(dist2mark < bestdist && thetaerror < thetabounds)
            bestdist = dist2mark;
            idx = ii;
        end
    end
end