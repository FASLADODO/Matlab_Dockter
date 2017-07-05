function XR = SortNorm(X)
%sorts data according to magnitude of norm in each row

    XD = NormRowWise(X);

    [~,idxs] = sort(XD);

    XR = X(idxs,:);
end