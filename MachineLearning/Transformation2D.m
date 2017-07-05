function T = Transformation2D(theta,shift,scale)

    T = [cos(theta), sin(theta), shift(1);
        -sin(theta),cos(theta),shift(2);
        0, 0, scale];
end