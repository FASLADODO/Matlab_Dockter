function [y, Jp] = scan (p)
% SCAN perform a range?and?bearing measure of a 2D point.
%
% In:
% p : point in sensor frame p = [p x ; p y]
% Out:
% y : measurement y = [range ; bearing]
% Jp: Jacobian wrt p

    px = p(1);
    py = p(2);
    d = sqrt(px^2+py^2);
    a = atan2(py,px);
    [a] = boundAngle(a);
    y = [d;a];
    
    if nargout > 1 % Jacobians requested
        Jp = [...
        px/sqrt(px^2+py^2) , py/sqrt(px^2+py^2)
        -py/(px^2*(py^2/px^2 + 1)), 1/(px*(py^2/px^2 + 1)) ];
    end
end