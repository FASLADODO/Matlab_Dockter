function [ground, moving, W, Z] = DyadSynthesisFunc(PP1,PP2,PP3,theta1,theta2,theta3,beta2,beta3)
    %This is merely a function version of the code used in Problem 2

    % setup
    deg2rad = (pi/180);

    %%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%
    % convert to radians because matlab
    theta1 = theta1 * deg2rad;
    theta2 = theta2 * deg2rad;
    theta3 = theta3 * deg2rad;
    beta2 = beta2 * deg2rad;
    beta3 = beta3 * deg2rad;

    % position differences (r*e^(i*theta) form )
    delta2 = abs(PP2 - PP1)*exp(i*angle(PP2 - PP1));
    delta3 = abs(PP3 - PP1)*exp(i*angle(PP3 - PP1));

    %angle differences
    alpha2 = theta2-theta1;
    alpha3 = theta3-theta1;


    % Cramers Rule Matrices
    denominator = [ exp(i*beta2) - 1, exp(i*alpha2) - 1;
                    exp(i*beta3) - 1, exp(i*alpha3) - 1 ];
    numerator1 = [ delta2, exp(i*alpha2) - 1;
                    delta3, exp(i*alpha3) - 1 ];
    numerator2 = [ exp(i*beta2) - 1, delta2;
                    exp(i*beta3) - 1, delta3 ];

    % Cramers rule for Ax=B using the 2 vector loop equations
    W = det(numerator1) / det(denominator);
    Z = det(numerator2) / det(denominator);

    %Ground Pivot
    ground = PP1 - Z - W;
    %moving pivot
    moving = PP1 - Z;


end
