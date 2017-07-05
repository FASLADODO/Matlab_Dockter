%%%%%%%%%%%%                 STU                 %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%% (Stewart-platform Totally Ultimate {3 DOF}) %%%%%%%%%%%%%%%%%%%%%%

%http://cdn.intechopen.com/pdfs-wm/34392.pdf

%radius of platforms
r_p     = 100; %radius top platform (mm)
r_base  = 75; %radius base platform(mm)

%Arm Lengths
A1 = 40;
A2 = 60;

%constants
deg2rad = (pi/180);
R0 = [0;0;0];

%input position/orientation
P_xyz     = [0;0;50]; %(x,y,z)
B_xyz     = [0;0;0]; %(x,y,z)


forward = -20:1:20;
backward = 20:-1:-20;
P_plot = [ forward,20*ones(1,length(forward)),backward,-20*ones(1,length(backward)) ;
           20*ones(1,length(forward)),backward,-20*ones(1,length(backward)),forward;
           50*ones(1,length(forward)),50*ones(1,length(forward)),50*ones(1,length(forward)),50*ones(1,length(forward))];
       
      
%get symbolic jacobian
[A] = JacobianStewart3(r_p,r_base,A1,A2);

for( ii = 1:length(P_plot) )
    
    P_xyz = P_plot(:,ii);
       
    %attachment point locations (base and platform)
    B(:,1) = [r_base;0;0];
    B(:,2) = [-r_base/2; (sqrt(3)/2)*r_base; 0];
    B(:,3) = [-r_base/2; -(sqrt(3)/2)*r_base; 0];

    P(:,1) = [r_p;0;0];
    P(:,2) = [-r_p/2; (sqrt(3)/2)*r_p; 0];
    P(:,3) = [-r_p/2; -(sqrt(3)/2)*r_p; 0];

    %unity rotation matrix
    R = [1,0,0;
        0,1,0;
        0,0,1];

    % Get link lengths/ servo angles
    L = [];
    theta = [];
    for i = 1:3
        L(:,i) = -B(:,i) + P_xyz + R*P(:,i);
        L_m(i) = norm(L(:,i));
        theta(i) = acos( (L_m(i)^2 + A1^2 - A2^2)/(2*A1*L_m(i)) ) - (pi/2);
    end

    %unit vectors
    U_p = ( B(:,1)+L(:,1)-P_xyz ) / norm(B(:,1)+L(:,1)-P_xyz);
    U_b = ( B(:,1) ) / norm(B(:,1));

    %circle plotter
    P_p = circle3D(r_p,P_xyz,R0,U_p,0);
    P_b = circle3D(r_base,B_xyz,R0,U_b,0);

    %plot it
    figure(1)
    plot3(P_p(1,:),P_p(2,:),P_p(3,:),'r','LineWidth',3);
    hold on
    plot3(P_b(1,:),P_b(2,:),P_b(3,:),'b','LineWidth',3);
    hold on
    PlotActuators(B,L);
    hold off
    
    axis([-150,150,-150,150,0,100])
    view(45,45)
    grid on
    title('STU - Stewart Platform 3')
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')


    PlotServos(theta,L_m,A1,A2);
    
    pause(0.3);

end

L
L_m
theta