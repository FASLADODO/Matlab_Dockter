function [X,Y,Z] = sphereObject(radius,center)
    
    %create it
    [X,Y,Z] = sphere(30); 
    
    %scale
    X = X.*radius;
    Y = Y.*radius;
    Z = Z.*radius;

    %Re-center
    X = X + center(1);
    Y = Y + center(2);
    Z = Z + center(3);
    
end