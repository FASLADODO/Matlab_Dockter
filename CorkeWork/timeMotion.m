%% Time stuff

%sin
[s,sdot,sdotdot] = tpoly(0,1,50);

plot(s)
hold on
plot(sdot)
hold on
plot(sdotdot)
hold off

smooth = mean(sdot)/max(sdot)

%% linear

[s,sdot,sdotdot] = lspb(0,1,50);

plot(s)
hold on
plot(sdot) %is a trapezoid
hold on
plot(sdotdot)
hold off

smooth = mean(sdot)/max(sdot)

%% multi dimensional

x= mtraj(@tpoly, [0,2], [1,-1], 50); 

plot(x)

%% segments

via = [4,1; 4,4; 5,2; 2,5];
q = mstraj(via, [2,1], [], [4,1], 0.05, 0 ); %create trajectory

plot(q)
hold on
scatter(via(:,1),via(:,2))
hold off

%% orientation interpolation

R0 = rotz(-1)*roty(-1);
R1 = rotz(1)*roty(1);

rpy0 = tr2rpy(R0); %roll pitch yaw
rpy1 = tr2rpy(R1);

rpy = mtraj(@tpoly, rpy0, rpy1, 50); %create trajectory

tranimate( rpy2tr(rpy));

%% rot+tran interpolation

T0 = transl(0.4,0.2,0) * trotx(pi);
T1 = transl(-0.4,-0.2,0.3) * troty(pi/2);

Ts = trinterp(T0, T1, [0:49]/49 );
about(Ts)

tranimate(Ts)

%% fancy interp
T0 = transl(0.4,0.2,0) * trotx(pi);
T1 = transl(-0.4,-0.2,0.3) * troty(pi/2);

Ts = ctraj(T0,T1,50);
P = transl(Ts);

plot(P)


%% initial velocity

initvel = 0.5;
finalvel = 10;

%sin
[s,sdot,sdotdot] = tpoly(0,1,50,initvel,finalvel);

% plot(s)
% hold on
 plot(sdot)
% hold on
% plot(sdotdot)
% hold off

smooth = mean(sdot)/max(sdot)









