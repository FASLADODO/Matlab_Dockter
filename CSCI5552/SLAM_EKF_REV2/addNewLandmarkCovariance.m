function P = addNewLandmarkCovariance(P,L,S)
    %P: current covariance
    %L: SLAM jacobians
    %S measurement noise
    %resize and fill in covariance matrix with new landmark
    
    %get current indices for all landmarks
    NN = size(P,1);
    nli = (NN - 3) / 2;
    newli = nli + 1; %new landmark index
    idh = 2 + (newli)*2;
    r = [1:3]; %robot indices
    l = [idh:idh+1]; %new landmark indices
    
    %resize covariance matrix
    Pnew = P;
    Pnew = [Pnew, zeros(NN,2)];
    Pnew = [Pnew; zeros(2,NN+2)];
    
    %covariance: Robot <-> new landmark
    G = P(r,r)*L.r';
    F = G';
    Pnew(r,l) = G;
    Pnew(l,r) = F;
    
    %covariance: new landmark <-> new landmark 
    C =  L.r * P(r,r) * L.r' + L.y * S * L.y';
    Pnew(l,l) = C;
    
    %landmark - old landmarks
    for ii = 1:nli
       %indices of existing landmarks
       idx = 2 + 2*ii;
       m = [idx:idx+1];
       
       %covariance: robot <-> existing landmarks
       Pri = P(r,m);
       Pri = L.r*Pri;
       
       %covariance: new landmark <-> existing landmarks
       Pnew(l,m) = Pri;
       Pnew(m,l) = Pri';
    end

    %store new covariance
    P = Pnew;
end