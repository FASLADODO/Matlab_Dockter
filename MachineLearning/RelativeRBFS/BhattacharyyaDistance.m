function D = BhattacharyyaDistance(p,q)
% In statistics, the Bhattacharyya distance measures the similarity of two discrete or continuous probability distributions.
%https://en.wikipedia.org/wiki/Bhattacharyya_distance

BC = sum(sqrt( p.*q ));
D = log(BC);

end