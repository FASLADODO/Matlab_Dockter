function [ T ] = DHFKRaven( jointvars, jointLengths )
    %Compute forward kinematics use DH method (For Raven)
    %Corke won't work since Lum uses craig convention
    
    
    %Get transformation matrices, craig convention
    T1 = Craig_DH(jointvars(1),0,0,0);
    T3 = Craig_DH(jointvars(2),0,0,-jointLengths.a13);
    T5 = Craig_DH(jointvars(3),0,0,jointLengths.a35);

    %get transformation matrix for given joint variables
    T = T1*T3*T5;

end

