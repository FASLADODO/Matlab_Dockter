%Corke transformations Chap 2

% se2 = [R,t;
%     0,1]; = 
%     [cos, sin, x; 
%     -sin,cos,y;
%     0,0,1];

%se2(X,Y,theta)

T1 = se2(1,2,30*pi/180) %transformation 1

T2 = se2(2,1,0) %transformation 2

T3 = T1*T2 %compose

T4 = T2*T1 %inverse compose

P = [3;2]; %random point

%P1 = h2e( inv(T1) * e2h(P) ); %h2e homogenous to euclidean
%or
P1 = homtrans(inv(T1), P) %more compact

P2 = homtrans(inv(T2), P) %relative to T2

trplot2(T1,'frame','1','color','b')
hold on
trplot2(T2,'frame','2','color','r')
hold on
trplot2(T3,'frame','3','color','g')
hold on
trplot2(T4,'frame','4','color','c')
hold on
plot_point(P,'*')

hold off

axis([0,5,0,5]);

%% 3D Rotations

Rx = rotx(pi/2) %rotation about x in 3d
Ry = roty(pi/2) %rotation about y in 3d
Rz = rotz(pi/2) %rotation about z in 3d

%tranimate(Rx) %the tightest shit

figure(1)
trplot(Rx,'frame','1','color','b')

figure(2)
trplot(Ry,'frame','2','color','r')

figure(3)
trplot(Rx*Ry,'frame','3','color','g')

figure(4)
trplot(Ry*Rx,'frame','4','color','c')

%% Euler angles

R = eul2r(0.1,0.2,0.3);
tranimate(R)

gamma = tr2eul(R)

%% Bryan angles

R = rpy2r(0.1,0.2,0.3);
[theta, v] = tr2angvec(R)

[v,lamba] = eig(R) %eigen values and vectors

%% quaternions q = scalar + vector = s <v1,v2,v3>

q = Quaternion( rpy2tr(0.1, 0.2, 0.3) )
q.norm
q.inv()

r = q.R

q.plot()

%% transformation

ty = transl(0,1,0) %translation matrix
tx = transl(1,0,0)

T = tx * trotx(pi/2) * ty

tranny = transl(T)

tranimate(T)



