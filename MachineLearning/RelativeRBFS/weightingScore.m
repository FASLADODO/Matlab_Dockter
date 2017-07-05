function metric = weightingScore(D)
%computes weighting of two parameter maximization

%See testweightingscore.m

% metric = min(D,[],2);
% metric = mean(D,2)*3 - abs(D(:,1) - D(:,2));
% metric = mean(D,2)*4 - abs(D(:,1) - D(:,2));
% metric = mean(D,2)*8 - abs(D(:,1) - D(:,2));
metric = mean(D,2);

end