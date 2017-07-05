function [Update, CurrentLandmarks] = localization2D(Prior, Landmarks, Measurements,matchthresh)
    %Prior = [X,Y,theta] current location guess in global
    %Landmarks = [X,Y; X, Y; ...] %Landmark locations in global
    %Measurements = [d, phi; d, phi; ...] %distance, bearing to visible
    %landmarks (local)
    %matchthresh is maximum distance for keeping points
    
    [NN,~] = size(Measurements);
    
    %Get currently visible landmarks in global space
    CurrentLandmarks = [];
    for ii = 1:NN
        %d,phi to X,Y
        temp = [Measurements(ii,1)*cos(Measurements(ii,2) + Prior(3)), Measurements(ii,1)*sin(Measurements(ii,2) + Prior(3)) ];
        %add in the prior to get global
        CurrentLandmarks(ii,:) = [temp(1) + Prior(1), temp(2) + Prior(2) ];
    end
    
    %Get all possible Landmark locations
    [kl,dist] = dsearchn(Landmarks,CurrentLandmarks);
    
    %create matchlist
    matchlist = [kl, [1:length(kl)]' ];
    
    %Weed out bad matches, match ID
    mid = matchlist(dist < matchthresh,:)
    
    %From here http://stackoverflow.com/questions/11687281/transformation-between-two-set-of-points
    Msub = [];
    proj = [];
    for jj = 1:length(mid)
        tempproj = Landmarks(mid(jj,1),:)';
        tempm = [ CurrentLandmarks(mid(jj,2),:) ,1,0;
                  CurrentLandmarks(mid(jj,2),2),-CurrentLandmarks(mid(jj,2),1),0, 1];
                
        Msub = [Msub; tempm];
        proj = [proj; tempproj];
    end
    trans = pinv(Msub)*proj;
    
    x_delta = trans(3); %M_13
    y_delta = trans(4); %M_23
    %acos(M11) asin(M12) -asin(M21) acos(M22)
    theta_delta = atan2(-trans(2),trans(1));
    %theta_delta = mean([acos(trans(1)), asin(trans(2)), -asin(trans(4)), acos(trans(5))] );
    %Update = Prior + [x_delta, y_delta, theta_delta];
    
    Update = [Prior(1)*trans(1)+Prior(2)*trans(2)+trans(3),-Prior(1)*trans(2)+Prior(2)*trans(1)+trans(4),Prior(3)+theta_delta];
    

end