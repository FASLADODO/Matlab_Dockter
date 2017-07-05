function [y, Jr, Jp] = observe(r, p)
% OBSERVE Transform a point P to robot frame and take a
% range?and?bearing measurement.
%
% In:
% r : robot frame r = [r x ; r y ; r alpha]
% p : point in global frame p = [p x ; p y]
% Out:
% y: range?and?bearing measurement
% Jr: Jacobian wrt r
% Jp: Jacobian wrt p

    [pr, PRr, PRp] = toFrame(r, p);
    [y, Jpr] = scan(pr);
    % The chain rule!
    Jr = Jpr * PRr;
    Jp = Jpr * PRp;

end