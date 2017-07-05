function [pf, Jf, Jp] = toFrame(F , p)
% TOFRAME transform point P from global frame to frame F
%
% In:
% F : reference frame F = [f x ; f y ; f alpha]
% p : point in global frame p = [p x ; p y]

% Out:
% pf: point in frame F
% Jf: Jacobian wrt F
% Jp: Jacobian wrt p

    t = F(1:2);
    a = F(3);
    R = [cos(a) -sin(a) ; sin(a) cos(a)];
    pf = R' * (p - t);
    
    if nargout > 1 % Jacobians requested
        px = p(1);
        py = p(2);
        x = t(1);
        y = t(2);
        Jf = [...
        [ -cos(a), -sin(a), cos(a)*(py - y) - sin(a)*(px - x)]
        [ sin(a), -cos(a), - cos(a)*(px - x) - sin(a)*(py - y)]];
        Jp = R';
    end
end
