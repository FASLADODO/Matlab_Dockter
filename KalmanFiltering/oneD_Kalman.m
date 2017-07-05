%test 1D kalman

%this is new

%one dimensional example
% https://en.wikipedia.org/wiki/Kalman_filter#Example_application

timesteps = 1000;
dt = 0.01;  %seconds

sigmaa = 0.9;
sigmaz = 0.9;

%initialize state and covariances
Xk1k1 = [0;0]; %start at zero
Pk1k1 = randn(2,2)*sigmaa;

%True position for simulation
Xv = 0.5; %m/s
Xtrue = [0; Xv];

xstore = [];

for tt = 1:timesteps
    %define our matrices
    Fk = [1, dt; 0, 1];
    G = [0.5*dt^2; dt];
    
    %compute noise
    w = G * normrnd(0,sigmaa);
    
    %predict
    Xkk1 = Fk*Xk1k1 + w;
    Pkk1 = F*Pk1k1*Fk';
    
    %update
    Xk = Xtrue + randn(2,1)*sigmaz; %simulated position update
    Hk = [1,0];
    Zk = Hk*Xk;
    R = sigmaz^2;
    
    %update equations
    yk = Zk - Hk*Xkk1;
    Sk = Hk*Pkk1*Hk' + R;
    Kk = Pkk1*Hk'*inv(Sk);
    Xkk = Xkk1 +Kk*yk;
    Pkk = (eye(2) - Kk*Hk)*Pkk1;
    
    %lets store
    xstore = [xstore; Xtrue(1), Xk(1), Xkk(1)];
    
    %next time step variables
    Xk1k1 = Xkk;
    Pk1k1 = Pkk;
    Xtrue = [1, dt; 0, 1] * Xtrue;
end

figure
plot(xstore(:,1),'r')
hold on
plot(xstore(:,2),'b')
hold on
plot(xstore(:,3),'c')
hold off
legend('true','observation','update')

