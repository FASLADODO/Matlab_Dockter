% Target region, area 5x5 cm
X = 5; % cm 
Y = 5; % cm 

% Make planar 2 link robot (symbolic version)
syms L1 L2 q1 q2
L(1) = Link( [ 0 0 L1 0]);
L(2) = Link( [ 0 0 L2 0]);
TwoLink = SerialLink(L, 'name', 'my2Link')
FKsym   = TwoLink.fkine ( [q1 q2] ) % NOTE: You can get symbolic FK & IK
IK = TwoLink.ikine_sym(2);          %  can also get symbolic IK for 2DOF
IK = IK{1}';                        %  just use one of the solutions    
Jsym    = TwoLink.jacob0( [q1 q2] ) % NOTE: can always get a symbolic Jacobian
Jsym    = Jsym(1:2, 1:2)            % 
det(Jsym)                           % and symbolic deteriminant

% You can even get symbolic singular values
disp('Get singular values (symbolic!):')
sigma = simplify( eig(Jsym.' * Jsym));
s1 = sigma(1) ;  % get the lowest singular value
pretty(sigma)


%% From lecture, optimal manipulability roughly at q2 = +/- pi/2; any q1
q1 = 0;
q2 = pi/2;
ss = eval(s1);   % substitute in in q1, q2
ss = vpa( ss,2)  % approximate to decimal

% plot a singular value over L1, L2 space
figure (100)
ezsurf( char(ss), [0 10, 0 10])
title('Singular Value of Jacobian')
xlabel('L1 [cm]')
ylabel('L2 [cm]')
zlabel('\sigma_1')


% To get your optimal solution, you can use this symbolic approach above 
% or a numerical approach similar to lecture (see below)

%% t = [0:.05:2]'; % time vector
L1 = 50;  % [cm] < ---- you need to pick these
L2 = 50;  % [cm]
L(1) = Link( [ 0 0 L1 0]);
L(2) = Link( [ 0 0 L2 0]);
TwoLink = SerialLink(L, 'name', 'my2Link')
q0 = [0 0];
qf = rand(1,2)*2;
%q = jtraj(q0, qf, t); % generate joint coordinate trajectory

figure(1)
TwoLink.plot(q0); hold on
% TwoLink.plot(qf);


% Jacobians from Corke:
% 
% % jacobn - Jacobian matrix in tool frame
% % jacob - dot Jacobian derivative
% % maniplty - manipulability
% % vellipse - display velocity ellipsoid
% % fellipse - display force ellipsoid
% % qmincon - null space motion to centre joints between limits

% figure(10)
% TwoLink.vellipse(qf); hold on
% TwoLink.vellipse(q0); hold on
%J0 = TwoLink.jacob0 (q0)
% 
% figure(11);clf; 
% for i = [0:.5:2*pi]
%     for j = [0:.2:2*pi]
% 
%         J = TwoLink.jacob0([i,j]);
%         %plot_ellipse( J(1:3)*J(1:3)'); hold on     
% 
%         % plot over joint space
%         plot3(i,j, TwoLink.maniplty([i j], 'yoshikawa'),'rx'); hold on
%         plot3(i,j, TwoLink.maniplty([i j], 'asada'),'gx'); hold on
%         plot3(i,j, sqrt(det(J(1:2,1:2)*J(1:2,1:2)')),'b.'); hold on
%         
%         % can also plot over cartesian space
%         T = TwoLink.fkine([i,j]);
%         x = T(1,end);
%         y = T(2,end);
%         % ...
%     end
% end
% grid on
% xlabel('q1')
% ylabel('q2')
% zlabel('Manipulability')
% legend('Yoshikawa', 'asada' ,'sqrt(det(J(1:2)*J(1:2)''');


%%
figure(12);clf; 
Jdet = [];
for i = [0:.5:2*pi]
    for j = [0:.2:2*pi]
        J = TwoLink.jacob0([i,j]);
        %plot_ellipse( J(1:3)*J(1:3)'); hold on 
        J = J(1:2, 1:2);
        S = svd(J);
        T = TwoLink.fkine([i,j]);
        x = T(1,end);
        y = T(2,end);
        subplot(2,1,1)
        plot3(i,j, min(S),        'rx'); hold on
        plot3(i,j, max(S),        'gx'); hold on
        plot3(i,j, min(S)./max(S),'b.'); hold on
        
        subplot(2,1,2)
        plot3(x,y, min(S),        'rx'); hold on
        plot3(x,y, max(S),        'gx'); hold on
        plot3(x,y, min(S)./max(S),'b.'); hold on
        
    end
end
subplot(2,1,1)
grid on
xlabel('q1')
ylabel('q2')
zlabel('Manipulability')

subplot(2,1,2)
grid on
xlabel('x')
ylabel('y')
zlabel('Manipulability')




%% You can plot over some set of joint angles (you must choose them)
% and hope it covers your region of interest, or you can do something
% more clever ...
% either way, be sure to plot/sketch/or overlay your target 5cm x 5cm
% square on your final plot
figure(20);clf; 
Jdet = [];
Npoints = 15
for q1 = [ ([0:Npoints]/Npoints - 0.5) ] * 2*X/L1  
    for q2 = [ ([0:Npoints]/Npoints - 0.5) ] * 2*Y/L2 +pi/2 
        
        J = TwoLink.jacob0([q1,q2]);
        %plot_ellipse( J(1:3)*J(1:3)'); hold on 
        J = J(1:2, 1:2);
        S = svd(J);
       
        subplot(2,1,1)
        %plot3(q1,q2, min(S),        'rx'); hold on; grid on
        %plot3(q1,q2, max(S),        'gx'); hold on
        plot3(q1,q2, min(S)./max(S),'b.'); hold on
        
        subplot(2,1,2)
        P = TwoLink.fkine([q1 q2]);
        tx = P(1,end);
        ty = P(2,end);
        %plot3(tx,ty, min(S),        'rx'); hold on; grid on
        %plot3(tx,ty, max(S),        'gx'); hold on
        plot3(tx,ty, min(S)./max(S),'b.'); hold on
        
    end
end
subplot(2,1,1)
grid on
xlabel('q1')
ylabel('q2')
zlabel('Manipulability')

subplot(2,1,2)
grid on
xlabel('x')
ylabel('y')
zlabel('Manipulability')
