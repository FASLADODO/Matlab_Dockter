function [ angle_out2, angle_out3 ] = GroundSpecification( delta2, delta3, R1, R2, R3, angle_in2, angle_in3, DyadType )
    % Solve for input/ coupler angles using ground specification methods

    % DyadType = 'Motion' or 'Timing'

    epsilon = 10^(-2);
    
    %  Change of variables for ground pivot specification, 
    % known beta, unknown alpha
    if(DyadType == 'Timing')
        D1 = R2*exp(i*angle_in3) - R3*exp(i*angle_in2);
        D2 = R3 - R1*exp(i*angle_in3);
        D3 = R1*exp(i*angle_in2) - R2;
    elseif (DyadType == 'Motion') % known alpha, unknown beta
        D1 = R3*exp(i*angle_in2) - R2*exp(i*angle_in3);
        D2 = -R3 + R1*exp(i*angle_in3);
        D3 = R2 - R1*exp(i*angle_in2);
    else
        error('Invalid dyad type!');
    end

    %magnitude of D vectors
    D1_mag = abs(D1);
    D2_mag = abs(D2);
    D3_mag = abs(D3);

    %angle of D vectors
    D1_ang = angle(D1);
    D2_ang = angle(D2);
    D3_ang = angle(D3);

    % solving for angle using D and tri values using equation given in handout 18, pt ii
    angle2_alt1 =  acos( (D3_mag^2 - D1_mag^2 - D2_mag^2)/(2*D1_mag*D2_mag))  + D1_ang - D2_ang;
    angle2_alt2 =  -acos( (D3_mag^2 - D1_mag^2 - D2_mag^2)/(2*D1_mag*D2_mag))  + D1_ang - D2_ang;

    %check for trivial (using acos range)
    angle2_alt1 = wrap2pi(angle2_alt1); %put in -pi to pi
    angle2_alt2 = wrap2pi(angle2_alt2); %put in -pi to pi
    % if solution 1 is trivial (== beta2) then use alternate (-acos)
    if(abs(angle2_alt1 - angle_in2) < epsilon)
        angle_out2 = angle2_alt2;
    else
        angle_out2 = angle2_alt1;
    end

    %solve for alpha3 using alpha2
    angle_out3 = angle( (-D1_mag*exp(i*D1_ang) - D2_mag*exp(i*D2_ang)*exp(i*angle_out2)) / (D3_mag*exp(i*D3_ang)));

end

