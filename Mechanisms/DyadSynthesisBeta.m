% Homework 3, Problem 3
% Analytical Dyad program for K and m point circles
% Rod Dockter, Oct. 2014

% setup
clear all
deg2rad = (pi/180);

%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%
% Specify all Beta2 values to try (angle between input positions 1 and 2) (deg)
beta2 = [25,300];

% Precision positions real+i*im (x,y)
PP1 = -279 + 152*i;
PP2 = -176 + 117*i;
PP3 = -105 + 36*i;

% Precision angles (deg)
theta1 = 22.87;
theta2 = 352.75;
theta3 = 287.61;



%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%
%Allocate beta3 values
beta3 = 0;
beta3rad = beta3 * deg2rad;
numPos = 361; %number of angles to try;

% convert to radians because matlab
theta1 = theta1 * deg2rad;
theta2 = theta2 * deg2rad;
theta3 = theta3 * deg2rad;
beta2rad = beta2 .* deg2rad;

% position differences (complex form )
delta2 = PP2-PP1;
delta3 = PP3-PP1;

%angle differences
alpha2 = theta2-theta1;
alpha3 = theta3-theta1;

% iterate through beta2 array
numBetas = length(beta2);
for bs = 1:numBetas
    %iterate through beta3 values
    for j = 1:numPos

        % Cramers Rule Matrices
        denominator = [ exp(i*beta2rad(bs)) - 1, exp(i*alpha2) - 1;
                        exp(i*beta3rad) - 1, exp(i*alpha3) - 1 ];
        numerator1 = [ delta2, exp(i*alpha2) - 1;
                        delta3, exp(i*alpha3) - 1 ];
        numerator2 = [ exp(i*beta2rad(bs)) - 1, delta2;
                        exp(i*beta3rad) - 1, delta3 ];

        % Cramers rule for Ax=B using the 2 vector loop equations
        W = det(numerator1) / det(denominator);
        Z = det(numerator2) / det(denominator);

        % pivot locations
        A0 = PP1 - Z - W;
        A1 = PP1 - Z;

        %store each ground and moving pivot into array (different rows for
        %beta2)
        groundPivots(bs,j) = A0;
        movingPivots(bs,j) = A1;

        %increment beta3
        beta3 = beta3 + 1;
        beta3rad = beta3 * deg2rad;
    end

end

% plot moving pivots, ground pivots and precision points
f1 = figure(1);

%plot all the beta circles
legendmatrix=cell(numBetas*2 + 3,1);
pidx = 1;
cmap = hsv(numBetas);
for bs = 1:numBetas
    plot(real(groundPivots(bs,:)),imag(groundPivots(bs,:)),'-','LineWidth',2,'Color',cmap(bs,:));
    hold on
    plot(real(movingPivots(bs,:)),imag(movingPivots(bs,:)),':','LineWidth',2,'Color',cmap(bs,:));
    hold on
    %add to legend for beta values
    legendmatrix{pidx} = strcat('Ground Pivots. B2 = ',num2str(beta2(bs)));
    pidx = pidx + 1;
    legendmatrix{pidx} = strcat('Moving Pivots. B2 = ',num2str(beta2(bs)));
    pidx = pidx + 1;
end
% plot precision points
plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold on
plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);

%format plots
axis equal
title('Center/Circle Point Circles.')
xlabel('X (real)')
ylabel('Y (imaginary)')
legendmatrix{pidx} = strcat('PP1');
legendmatrix{pidx+1} = strcat('PP2');
legendmatrix{pidx+2} = strcat('PP3');
legend(legendmatrix,'Location','BestOutside')


%%
% Additional Plotting for Problem 3, Part D

%values found for part C
A0 = -251.937948 + -42.893825 *i 
A1 = -223.558398 + -10.210592 *i 
B0 = -298.165343 + 185.315481 *i 
B1 = -324.520955 + 51.833125 *i

% plot moving and ground point and reast of mechanism
plot(real(A0),imag(A0),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(A1),imag(A1),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B0),imag(B0),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(B1),imag(B1),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot([real(A0), real(A1)],[imag(A0), imag(A1)],'-','LineWidth',4,'Color',[0 0 0]);
hold on
plot([real(B0), real(B1)],[imag(B0), imag(B1)],'-','LineWidth',4,'Color',[0 0 0]);
hold on
plot([real(B1), real(PP1)],[imag(B1), imag(PP1)],'-','LineWidth',4,'Color',[0 0 0]);
hold on
plot([real(A1), real(PP1)],[imag(A1), imag(PP1)],'-','LineWidth',4,'Color',[0 0 0]);
hold on
plot([real(A1), real(B1)],[imag(A1), imag(B1)],'-','LineWidth',4,'Color',[0 0 0]);

hold off


%%
% GUI for clicking on one circle and getting point on other

gg = 1;
while 1
    [x,y] = ginput(1);
    
    %delete last plotted point
    if gg > 1
       delete(mplot); 
       delete(gplot); 
    end

    % sorts through click point to find if its nearest a ground or moving
    % circle
    tmpg = sqrt((real(groundPivots(1,:))-x).^2 + (imag(groundPivots(1,:))-y).^2);
    [sg idxg] = min(tmpg);
    tmpm = sqrt((real(movingPivots(1,:))-x).^2 + (imag(movingPivots(1,:))-y).^2);
    [sm idxm] = min(tmpm);
    

    % based on closest circle specifies point for both pivots
    if(sg < sm)
        closestm = movingPivots(1,idxg) ;
        closestg = groundPivots(1,idxg) ;
    else
        closestm = movingPivots(1,idxm) ;
        closestg = groundPivots(1,idxm) ;
    end
    
    % plot moving and ground point
    mplot = plot(real(closestg),imag(closestg),'o','MarkerSize',3,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
    hold on
    gplot = plot(real(closestm),imag(closestm),'o','MarkerSize',3,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
    %increment indices
    gg = gg + 1;
    
    %prints
    fprintf(' Ground Pivot = %f + %f *i \n',real(closestg), imag(closestg));
    fprintf(' Moving Pivot = %f + %f *i \n',real(closestm), imag(closestm));
    
    %if mouse clicked go again, if keyboard press, exit
    w = waitforbuttonpress;
    if w == 0
        continue;
    else
        break;
    end
end

hold off



