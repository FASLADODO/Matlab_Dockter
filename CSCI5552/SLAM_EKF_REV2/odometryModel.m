function [odomNext,Jr,Jn] = odometryModel(odomLast,u,q)
%IN:
% odomLast = [x,y,theta]
% u: control signal u = [d x ; d alpha] = [time step * velocity; time step * anglediff
% q = noise parameter
% Out:
% ro: updated robot pose
% Jr: Jacobian d(ro) / d(r)
% Jn: Jacobian d(ro) / d(n

    %figure out input
    a = odomLast(3);
    dx = u(1) + q(1)*randn(1,1);
    da = u(2) + q(2)*randn(1,1);
    
    ao = a + da;
    %bound angle
    ao = boundAngle(ao);
    
    dp = [dx;0];
    
    %get jacobians and update position in frame
    [to, T0r, T0dt] = fromFrame(odomLast, dp);
    A0a = 1;
    A0da = 1;
    Jr = [T0r ; 0 0 A0a];
    Jn = [T0dt(:,1) zeros(2,1) ; 0 A0da];
    
    odomNext = [to;ao];
end