function scan = generate_scan(x, y, phi, sigma_r)
% THIS FUNCTION PLOTS A SET OF OBSTACLES (SPECIFIED BY THE USER) 
% AND PRODUCES THE REAL AND MEASURED LASER SCAN FROM A LOCATION (X, Y, PHI)
%
% Output Variables
% scan.r:       True distances
% scan.rm:      Measured distances
% scan.theta:   True Bearing-angles
%
% Input Variables
% x, y:         Position of Robot in Global Frame (m)
% phi:          Bearing of Robot in Global Frame (rad)
% sigma_r:      Range msmt. standard deviation (m)


% Flag whether to plot figures (1) or not (0)
% depending on your system, this might make a big difference in speed
visualize = 0;

%---------------------------------------------------------------------
% LASER SCANNER INTRINSIC PARAMETERS
%---------------------------------------------------------------------

% Laser scanner range
RANGE = 8; % meters

% Laser beam angular seperation
laser_delta_theta=pi/360; % radians

% Standard Deviation of the noise in the laser distance measurements
std_laser= sigma_r; % meters

% start and end of the field of view of the laser scanner
phiL_start = -pi/2; % radians
phiL_end = pi/2;    %radians

%---------------------------------------------------------------------



%---------------------------------------------------------------------
% LASER SCANNER EXTRINSIC PARAMETERS (from input)
%---------------------------------------------------------------------

% POSITION & ORIENTATION
xL0 = x; % meters
yL0 = y; % meters
phiL0 = phi; % radians



% HOMOGENEOUS TRANSFORMATION MATRIX

% rotational matrix
R= [cos(phiL0) -sin(phiL0);
    sin(phiL0) cos(phiL0)];

p = [xL0;
     yL0];

Tr= [R          p;
    zeros(1,2) 1];

%---------------------------------------------------------------------
% OBSTACLE EXAMPLES *** YOU CAN SET YOUR OWN OBSTACLES ***
% Change N according to the number of obstacles you want to include
% Change M(i) to reflect the number of sides of obstacle i
%---------------------------------------------------------------------

% Number of POLYGONAL OBSTACLES 
N=2;

% M is the vector that determines the number of sides of each 
% polygonal obstacle

% OBSTACLE 1
M(1)=3;
x(1,1) = 2; y(1,1) = 2;
x(1,2) = 3; y(1,2) = 3;
x(1,3) = 4; y(1,3) = 2;


% OBSTACLE 2
M(2)=4;
x(2,1) = -2; y(2,1) = 4;
x(2,2) =  1; y(2,2) = 4;
x(2,3) =  1; y(2,3) = 6;
x(2,4) = -2; y(2,4) = 6;

%---------------------------------------------------------------------



%---------------------------------------------------------------------
% PLOT OBSTACLES

if visualize
    figure(1); hold on;

    for i=1:N,
        plot([x(i,1:M(i)) x(i,1)], [y(i,1:M(i)) y(i,1)],'LineWidth',2);
    end
    xlabel('x (m)');
    ylabel('y (m)');

    axis equal;



    figure(2); hold on;

    for i=1:N,
        plot([x(i,1:M(i)) x(i,1)], [y(i,1:M(i)) y(i,1)],'LineWidth',2);
    end
    xlabel('x (m)');
    ylabel('y (m)');

    axis equal;
end

%---------------------------------------------------------------------





%---------------------------------------------------------------------
% Parameters of the Line Equations a*x + b*y = c for the 
% obstacle line segments:
%---------------------------------------------------------------------

for i=1:N,
   for j=1:M(i),
      if j<M(i)
         [a(i,j),b(i,j),c(i,j)]=line_equation(x(i,j),y(i,j),x(i,j+1),y(i,j+1));
      else
         [a(i,j),b(i,j),c(i,j)]=line_equation(x(i,M(i)),y(i,M(i)),x(i,1),y(i,1));
      end
   end
end

%---------------------------------------------------------------------



%---------------------------------------------------------------------
% LASER BEAMS STARTING, ENDING POINTS AND LINE PARAMETERS
%---------------------------------------------------------------------

% Number of laser beams
NL = floor(abs(phiL_end-phiL_start)/laser_delta_theta)+1;


for k=1:NL,
    % Starting point of laser beam i in local coordinates
    L_p=[0;
         0;
         1];

    % Starting point of laser beam i in global coordinates 
    G_p = Tr*L_p;
    xL(k,1)= G_p(1,1);
    yL(k,1)= G_p(2,1);   

    % Ending point of laser beam i in local coordinates
    L_p=[RANGE * cos(phiL_start + (k-1)*laser_delta_theta);
         RANGE * sin(phiL_start + (k-1)*laser_delta_theta);
         1];

    % Ending point of laser beam i in global coordinates 
    G_p = Tr*L_p;
    xL(k,2)= G_p(1,1);
    yL(k,2)= G_p(2,1);   

    % LINE PARAMETERS
    [aL(k),bL(k),cL(k)]=line_equation(xL(k,1),yL(k,1),xL(k,2),yL(k,2));
end



% Add laser lines to plot
if visualize
    figure(1);
    hold on;

    for k=1:NL,
        plot(xL(k,:),yL(k,:),'r:','LineWidth',0.1);
    end

    figure(1); hold off;
end

%---------------------------------------------------------------------
% LASER BEAMS INTERSECTION WITH OBSTACLE SIDES
%---------------------------------------------------------------------


% for every laser beam
for k=1:NL,
    % initially assume that no point of intersection will be found
    count=0;
    % Vectors used to store the length and coordinates of all points intersected 
    % by a laser beam
    store_length_OT=0;
    storeTx=0;
    storeTy=0;

    % for every obstacle
    for i=1:N,
        % for every obstacle side
        for j=1:M(i),
            D=[aL(k)  bL(k);
                a(i,j) b(i,j)];

            % check if the a line segment AB and a laser beam OP intersect
            if (abs(det(D))<0.01)
                % the two lines are parallel - no intersection point
            else
                % the two lines intersect at the following point
                % T = inv(D)*[cL(k); c(i,j)];
                % T: point of intersection
                T = [( b(i,j)*cL(k) - bL(k)*c(i,j) )/det(D);
                    ( - a(i,j)*cL(k) + aL(k) * c(i,j) )/det(D)];                

                % O: starting point of the laser beam
                O= [xL(k,1);
                    yL(k,1)];

                % distance to the point of intersection from the laser source
                OT = T - O;
                sq_length_OT = OT'*OT;

                % Check if out of range
                if (sq_length_OT > RANGE^2)
                    % Point T is out of range
                else
                    % Point T is within range
                    %check if this point is within the obstacle line segment
                    % A: start of obstacle line segment
                    A=[x(i,j);
                       y(i,j)];

                    % B: end of obstacle line segment
                    if (j == M(i))
                        B=[x(i,1);
                            y(i,1)];
                    else
                        B=[x(i,j+1);
                            y(i,j+1)];
                    end

                    AT=T-A;
                    length_AT = sqrt(AT'*AT);

                    TB=B-T;
                    length_TB = sqrt(TB'*TB);

                    AB=B-A;
                    length_AB = sqrt(AB'*AB);

                    if (abs(length_AT + length_TB - length_AB)>0.01)
                        % this point intersects the extension of segment AB
                    else
                        % point within segment AB
                        % check if this point is within the beam line segment
                        % P: ending point of the laser beam
                        P= [xL(k,2);
                            yL(k,2)];

                        % T point of intersection as before
                        length_OT = sqrt(sq_length_OT);

                        TP= P-T;
                        length_TP = sqrt(TP'*TP);

                        length_OP=RANGE;

                        if (abs(length_OT+length_TP-length_OP)>0.01)
                            % this point intersects with an extension of the segment OP
                        else
                            % point within the beam OP
                            count=count+1;
                            store_length_OT(count)=length_OT;
                            storeTx(count)=T(1,1);
                            storeTy(count)=T(2,1);
                        end
                    end
                end
            end
        end % for j
    end % for i

    % check to find if any points were found
    if (count>0)
        % if more than one points of intersection found select the closest
        [value index]=min(store_length_OT);
        T=[storeTx(index);
           storeTy(index)];

       % d is the distance measured by the laser scanner
       d(k)= value;
    else
        % if no point of intersection was found return defaul valus of 0
        T=[xL0;
           yL0];

       % d is the distance measured by the laser scanner
       d(k)=0;
    end

    % plot Intersection point
    if visualize
        figure(2); hold on; plot(T(1,1),T(2,1),'kx','LineWidth',2);
    end

end % for k


%---------------------------------------------------------------------
% LASER BEAMS WITH ADDITIVE GAUSSIAN NOISE
%---------------------------------------------------------------------

% measured distance dm
dm=zeros(1,NL);

for k=1:NL,
    if (d(k)==0)
        % no return add no noise
    else
        dm(k)= d(k) + randn * std_laser;
    end

    % Intersection point of laser beam i in local coordinates
    L_p=[dm(k) * cos(phiL_start + (k-1)*laser_delta_theta);
         dm(k) * sin(phiL_start + (k-1)*laser_delta_theta);
         1];

    % Intersection point of laser beam i in global coordinates 
    G_p = Tr*L_p;
    xLm(k,2)= G_p(1,1);
    yLm(k,2)= G_p(2,1);


    % plot measured intersection point
    if visualize
        figure(2); hold on;
        plot(xLm(k,2),yLm(k,2),'r+');
    end

end



if visualize
    figure(2); hold off;

    % Plot true and measured laser scan, distance vs. bearing
    figure; hold on;
    plot( [phiL_start:laser_delta_theta:phiL_end] , d  ,'kx');
    plot( [phiL_start:laser_delta_theta:phiL_end] , dm ,'r+');
    xlabel('angle (rad)');
    ylabel('distance (m)');
    legend('Real Distance','Measured Distance');

    hold off;
end


% fill in return value
scan.r  = d;  % true distances
scan.rm = dm; % measured distances
scan.theta = phiL_start : laser_delta_theta : phiL_end; % measured bearigns (noise-free, actually)