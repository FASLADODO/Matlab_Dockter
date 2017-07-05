% Analytical 4 position motion generation (Burmeister curve)
% Rod Dockter, Dec. 2014
% using handout # 23 Dr. Sezen notes

% setup
clear all
deg2rad = (pi/180);

%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%

% Precision positions real+i*im (x,y)
PP1 = -309 + 143*i;
PP2 = -186 + 69*i;
PP3 = -78 +- 12*i;
PP4 = -2 - 124*i;

% Precision angles (deg)
theta1 = 32.0;
theta2 = 5.8;
theta3 = 341.9;
theta4 = 326.0;


%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%

% convert to radians because matlab
theta1 = theta1 * deg2rad;
theta2 = theta2 * deg2rad;
theta3 = theta3 * deg2rad;
theta4 = theta4 * deg2rad;

% position differences (complex form )
delta2 = PP2-PP1;
delta3 = PP3-PP1;
delta4 = PP4-PP1;

%angle differences
alpha2 = theta2-theta1;
alpha3 = theta3-theta1;
alpha4 = theta4-theta1;

%Allocate beta2 values
beta2 = 0;
beta2rad = beta2 * deg2rad;
numPos = 361;

%iterate through beta3 values
for jj = 1:numPos

    %D variables
    D3 = delta2*(exp(i*alpha4) - 1) - delta4*(exp(i*alpha2) - 1); % doesnt change
    D4 = delta3*(exp(i*alpha2) - 1) - delta2*(exp(i*alpha3) - 1); % doesnt change
    D = (delta3 - delta2)*exp(i*alpha4) + (delta4 - delta3)*exp(i*alpha2) + (delta2 - delta4)*exp(i*alpha3) + (delta4*(exp(i*alpha3) - 1 ) - delta3*(exp(i*alpha4) - 1) ) *exp(i*beta2rad);

    % angles and mags
    argD3 = angle(D3);
    argD4 = angle(D4);
    argD = angle(D);
    magD3 = abs(D3);
    magD4 = abs(D4);
    magD = abs(D);
    
    % get missing betas 
    beta4 = [ acos( (magD3^2 - magD^2 -magD4^2) / (2*magD*magD4) ) + argD - argD4 , -acos( (magD3^2 - magD^2 -magD4^2) / (2*magD*magD4) ) + argD - argD4 ];
    beta3 = [ angle( (-magD*exp(i*argD) - magD4*exp(i*argD4)*exp(i*beta4(1)) ) / (magD3*exp(i*argD3) ) ), angle( (-magD*exp(i*argD) - magD4*exp(i*argD4)*exp(i*beta4(2)) ) / (magD3*exp(i*argD3) ) ) ];

    for kk = 1:2
        % cramers matrix
        A = [ exp(i*beta2rad) - 1, exp(i*alpha2) - 1;
            exp(i*beta3(kk)) - 1, exp(i*alpha3) - 1;
            exp(i*beta4(kk)) - 1, exp(i*alpha4) - 1 ];
        B = [ delta2;
            delta3;
            delta4];

        % Cramers rule for Ax=B using the 2 vector loop equations
        XX = A\B;
        W = XX(1);
        Z = XX(2);
        
        % pivot locations
        A0 = PP1 - Z - W;
        A1 = PP1 - Z;

        %store each ground and moving pivot into array (different rows for
        %beta2)
        groundPivots(kk,jj) = A0;
        movingPivots(kk,jj) = A1;
        
    end

    %increment beta3
    beta2 = beta2 + 1;
    beta2rad = beta2 * deg2rad;
end

figure(1)
plot(real(groundPivots(1,:)),imag(groundPivots(1,:)),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(groundPivots(2,:)),imag(groundPivots(2,:)),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(movingPivots(1,:)),imag(movingPivots(1,:)),'x','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold on
plot(real(movingPivots(2,:)),imag(movingPivots(2,:)),'x','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold off
