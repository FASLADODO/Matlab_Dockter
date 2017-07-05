function PIN = MinBoundEllipse_Test(A,C,P)
%Returns percent of points within the ellipse
%A: covariance matrix from MivVolEllipse()
%C: centroid from  MivVolEllipse()
%P: Data Points to check
%PIN: percent inside

    for ii = 1:length(P)
        inside(ii) = (P(ii,:)' - C)'*A*(P(ii,:)' - C) <= 1;
    end

    PIN = sum(inside)/length(inside);

end