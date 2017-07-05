function [theta] = boundAngle(theta)
    if theta > pi
        theta = theta - 2*pi;
    end
    if theta < -pi
        theta = theta + 2*pi;
    end
end