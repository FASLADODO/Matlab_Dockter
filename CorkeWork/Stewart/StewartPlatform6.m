%%%%%%%%%%%%                 STU                 %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%% (Stewart-platform Totally Ultimate) %%%%%%%%%%%%%%%%%%%%%%

%http://cdn.intechopen.com/pdfs-wm/34392.pdf

offset = [0;0;50];
deg2rad = (pi/180);

%input position/orientation
P     = [0;50;0]; %(x,y,z)
P_dot = [0.2;-0.3;10];
R_xyz = [0;0.1;0]; %(alpha,beta,gamma)
R_dot = [0.01;0.02;0.03];

P = P+offset;

%radius
r_p     = 50; %radius top platform (mm)
r_base  = 100; %radius base platform(mm)
theta_p = [90*deg2rad;30*deg2rad]; %moving platform spacing (even=60)
theta_b = [90*deg2rad;30*deg2rad]; %base platform spacing (even=60)

%platform seperations
lambda = [];
for i = [1,3,5]
    lambda(i) = (i*pi/3) - (theta_p(1)/2);
end
for i = [2,4,6]
    lambda(i) = lambda(i-1) + theta_p(2);
end

%base seperations
v = [];
for i = [1,3,5]
    v(i) = (i*pi/3) - (theta_b(1)/2);
end
for i = [2,4,6]
    v(i) = v(i-1) + theta_b(2);
end

%attachment point locations
GT = [];
B = [];
for i = 1:6
    GT(:,i) = [r_p*cos(lambda(i)); r_p*sin(lambda(i)); 0];
    B(:,i)  = [  r_base*cos(v(i));   r_base*sin(v(i)); 0];
end

%rotation matrix
R_BT = [cos(R_xyz(2))*cos(R_xyz(3)), cos(R_xyz(3))*sin(R_xyz(1))*sin(R_xyz(2)) - cos(R_xyz(1))*sin(R_xyz(3)), sin(R_xyz(1))*sin(R_xyz(3)) + cos(R_xyz(1))*cos(R_xyz(3))*sin(R_xyz(2)) ;
        cos(R_xyz(2))*sin(R_xyz(3)), cos(R_xyz(1))*cos(R_xyz(3)) + sin(R_xyz(1))*sin(R_xyz(2))*sin(R_xyz(3)), cos(R_xyz(1))*sin(R_xyz(2))*sin(R_xyz(3)) - cos(R_xyz(3))*sin(R_xyz(1)) ;
                     -sin(R_xyz(2)),                                             cos(R_xyz(2))*sin(R_xyz(1)),                                             cos(R_xyz(1))*cos(R_xyz(2))];
                 
%get link lengths             
L_v = []; %vector
L_m = []; %magnitude
for i = 1:6
    L_v(:,i) = R_BT*GT(:,i) + P - B(:,i);
    L_m(i)   = norm(L_v(:,i));
    %L_m(i) = sqrt( (P(1) - B(1,i) + GT(1,i)*R_BT(1,1) + GT(2,i)*R_BT(1,2))^2 + (P(2) - B(2,i) + GT(1,i)*R_BT(2,1) + GT(2,i)*R_BT(2,2))^2 + (P(3) + GT(1,i)*R_BT(3,1) + GT(2,i)*R_BT(3,2))^2 ); %more optimal way
end

%unit vectors
U_p = ( B(:,1)+L_v(:,1)-P ) / norm(B(:,1)+L_v(:,1)-P);
U_b = ( B(:,1) ) / norm(B(:,1));

P_p = circle3D(r_p,P,R_xyz,U_p,0);
P_b = circle3D(r_base,[0;0;0],[0;0;0],U_b,0);

plot3(P_p(1,:),P_p(2,:),P_p(3,:),'r','LineWidth',3);
hold on
plot3(P_b(1,:),P_b(2,:),P_b(3,:),'b','LineWidth',3);
hold on
PlotActuators(B,L_v);
hold off

grid on
axis([-100,100,-100,100,0,100])
title('STU - Stewart Platform')
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')

L_m
%% Jacobians

% u_i = L_v_i/L_m_i
for ii = 1:6
   u(:,ii) = L_v(:,ii)./L_m(ii);
end


J_IB = zeros(6,6);

for ii = 1:6
   J_IB(ii,1) = u(1,ii); 
   J_IB(ii,2) = u(2,ii);
   J_IB(ii,3) = u(3,ii);
   
   prod = cross( R_BT*GT(:,ii),u(:,ii) );
   J_IB(ii,4:6) = prod';
end

J_IIB = [1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        0,0,0,cos(R_xyz(2)),0,0;
        0,0,0,0,1,-sin(R_xyz(1));
        0,0,0,-sin(R_xyz(2)),0,cos(R_xyz(1))];
    
X_dot = [P_dot;R_dot]; %cartesian Velocities

J = J_IB*J_IIB;

L_dot = J\X_dot;





