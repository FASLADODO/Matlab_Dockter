function mask = checkInsideHull(DT,Data)
%check if data is inside a convex hull or not
%DT: deluanay tirangulation object
%Data: N-D data matrix
%mask: 1 if inside, 0 if outside

simplexIndex = pointLocation(DT,Data);
mask = ~isnan(simplexIndex);

end