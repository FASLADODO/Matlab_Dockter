function [ T ] = DHFKRhino( jointvars )
    %Compute forward kinematics use DH method (For Rhino)
    
    %joint lengths
    d1 = 195;
    a2 = 170; 
    a3 = 170;
    a4 = 1; 
    d5 = 125;

    %offset zero pose, so it makes robot straight up in air, matches screw
    %formulation
    theta_offset = [0,-pi/2,0,-pi/2,pi]; %normally would be in DH table

    %theta,d,a,alpha (DH convention)
    L(1) = Link([0 d1 0 -pi/2]); %link 1
    L(2) = Link([0 0 a2 0]); %link 2
    L(3) = Link([0 0 a3 0]); %link 3
    L(4) = Link([0 0 a4 -pi/2]); %link 4
    L(5) = Link([0 d5 0 0]); %link 5

    %create robot object
    Rhino = SerialLink(L, 'name', 'Rhino');

    %get transformation matrix for given joint variables
    T = Rhino.fkine(jointvars+theta_offset);

end

