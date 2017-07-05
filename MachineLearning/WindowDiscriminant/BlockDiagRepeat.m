function M = BlockDiagRepeat(X,n)
% eg:
%     n = length(data);
%     X = 1./cov(data);
    Ac = repmat({X},n,1);
    M = blkdiag(Ac{:});

end