function L = UnitVec(x)
    %return the unit vector for any N dimensional vector x
    %x = [X,Y,Z] row vector
    L = x ./ norm(x);
end