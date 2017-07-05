function [match1, match2] = ICP_match(R1_p_R2_init, R1_phi_R2_init, scan1, scan2)
%%Rod Dockter
%%Assignment 2, problem 4
%%Program for matching points from scan 1 to scan 2 at each step
%%I used a pair of loops to go through scan1.rm and scan2.rm and search for
%%the non zero entries in the measured distances which correspond to actual
%%recorded objects. A second array, scan_pnt() records in order which
%%indices of scan.rm returned non-zero results.
%%Also in these loops I record the x and y cartesian coordinates
%%corresponding to each non-zero r and theta.

%%first to loops for identifying non-zero entries in scan1
%%Note pre-allocating is not an option so these loops may run slow
i = 1;
for j = 1:length(scan1.rm)
    if scan1.rm(j)>0
        scan1pnts(i) = j;
        [xscan1(i), yscan1(i)] = pol2cart(scan1.theta(j), scan1.rm(j));
        i=i+1;
    end
end
%%And non-zero entries in scan2
i = 1;
for j = 1:length(scan2.rm)
    if scan2.rm(j)>0
        scan2pnts(i) = j;
        [xscan2(i), yscan2(i)] = pol2cart(scan2.theta(j), scan2.rm(j));
        i=i+1;
    end
end

%%This second portion takes the arrays of the nonzero pnts of scan and runs
%%through a triple nested loop to find the match in scan2 that corresponds
%%to a minimization of the Mahalanobis distance.
%%To do this I used a first loop to run through the smaller of the two
%%non-zero scans in order to populate the match1 and match2 arrays
%%properly.
%%The next two loops run through the non-zero scan1 and scan2 points
%%respectively and if neither current scan's point has been flagged as
%%used, then it tries that combination of scan1 and scan2 points to check
%%if those two points are a succesful match using the normalized distance
%%metric.
%%This way any possible match is tried so the matches are not biased to the
%%first match.
%%When two scan points do get used those points are flagged in the
%%scan_matchpnts arrays as a 0 so they do not get used again, making sure
%%that the winner-takes-all strategy is followed.

%%Pre-allocating arrays for flagging used scanned pnts.
scan1_matchpnts = ones(length(scan1pnts),1);
scan2_matchpnts = ones(length(scan2pnts),1);

%%Determining size of matched arrays.
smallerscan = min(length(scan1pnts),length(scan2pnts));

%%Preformulating the rotation matrix
C_phi = [cos(R1_phi_R2_init),-sin(R1_phi_R2_init);sin(R1_phi_R2_init),cos(R1_phi_R2_init)];
fillm = 1;
%%looping trough potential matches
for m = 1:smallerscan
    previousmin = 10000;
    for n = 1:length(scan1pnts)
        for o = 1:length(scan2pnts)
            if scan1_matchpnts(n) == 1 && scan2_matchpnts(o) == 1
                %%Computing the normalized distances from equation given
                R1_z_li=[xscan1(n); yscan1(n)];
                R2_z_lj=[xscan2(o); yscan2(o)];
                newmin = norm(R1_z_li - (R1_p_R2_init + C_phi*R2_z_lj));
                if newmin < previousmin
                    previousmin = newmin;
                    scan1min = n;
                    scan2min = o;
                end
            end
        end
    end
    scan1_matchpnts(scan1min) = 0;
    scan2_matchpnts(scan2min) = 0;
    minarray(m) = previousmin;
    %%records indices of matches in scan1 and scan2
    %%Last I check to see if the distance metric for matching is diverging
    %%aka if all the good matches are gone
    %%if this is the case I do not record the new one as a match
    if m == 1
        match1(fillm) = scan1pnts(scan1min);
        match2(fillm) = scan2pnts(scan2min);
        fillm = fillm + 1;
    end
    %%Case when the change in the distance metric meets our threshold
    if m > 1
       if abs((minarray(m) - minarray(m-1))/2) < 0.0005
           match1(fillm) = scan1pnts(scan1min);
           match2(fillm) = scan2pnts(scan2min);
           fillm = fillm + 1;
       else
           fillm = fillm; %%Threshold not met.
       end
    end
end


end


