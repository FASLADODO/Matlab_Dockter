%Discriminant Least Squares test, based on talk 3-16-2015
%see pdf notes

clear all
close all
clc

%|xdot    | = |0     1|*|x   | + U
%|xdotdot |   |-a1 -a2| |xdot|

%Different Parameter vectors
paramOptions = 2;
phi1 = [100;30;1];
phi2 = [100;30;1];
  phi2 = [15;200;1];

type = {'params1','params2'};

%%%%%%%%%%%%%%%%%%%%%%start with first phi%%%%%%%%%%%%%%%%%%%%%%
phi = phi1;

A = 100;
tend = 10.25;
T = 0.001; % sampling period is fronm 1KHz
t = 0:0.001:tend;

%input is linear increase with time
input.time = t;
input.signals.values = A*t;
 
input.time = [input.time]';
input.signals.values = [input.signals.values]';
input.signals.dimensions = 1;

% SIMULINK
fprintf('Running simulation 1...')

sim('SimpleModel1.slx');
fprintf(' DONE\n')

% no noise, get output stuff
u1 = input_out.Data(:,1);
x1 = state.Data(:,3);
xdot1 = state.Data(:,2);
xdotdot1 = state.Data(:,1); 
state1 = state.Data;

%add in some noise
x1 = x1 + (rand(size(x1))-.5)*0.01*max(x1); 
xdot1 = xdot1 + (rand(size(xdot1))-.5)*0.01*max(xdot1);
xdotdot1 = xdotdot1 + (rand(size(xdotdot1))-.5)*0.01*max(xdotdot1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%now with second phi%%%%%%%%%%%%%%%%%%%%%%
phi = phi2;

% SIMULINK
fprintf('Running simulation 2...')

sim('SimpleModel1.slx');
fprintf(' DONE\n')

% no noise, get output stuff
u2 = input_out.Data(:,1); 
x2 = state.Data(:,3); 
xdot2 = state.Data(:,2); 
xdotdot2 = state.Data(:,1); 
state2 = state.Data;

%add in some noise
x2 = x2 + (rand(size(x2))-.5)*0.01*max(x2); 
xdot2 = xdot2 + (rand(size(xdot2))-.5)*0.01*max(xdot2);
xdotdot2 = xdotdot2 + (rand(size(xdotdot2))-.5)*0.01*max(xdotdot2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot both%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(t,state1(:,3),'r');
hold on
plot(t,state1(:,2),'g');
hold on
plot(t,state1(:,1),'b');
hold off
legend('x','xdot','xdotdot')
title('linear model, params 1')
xlabel('time')
ylabel('states')

figure(2)
plot(t,state2(:,3),'r');
hold on
plot(t,state2(:,2),'g');
hold on
plot(t,state2(:,1),'b');
hold off
legend('x','xdot','xdotdot')
title('linear model, params 2')
xlabel('time')
ylabel('states')


%% now use total least squares trick to back compute input1 and input2

%create data matrices
D1 = [x1,xdot1,xdotdot1];
D2 = [x2,xdot2,xdotdot2];

%TLS (should return original u1 and u2
u1_tls = D1*phi1;
u2_tls = D2*phi2;

%plot to check
figure(3)
plot(t,u1_tls,'r');
hold on
plot(t,u2_tls,'g');
hold off
legend('u1 tls','u2 tls')
title('Input, backcomputed from TLS')
xlabel('time')
ylabel('input val')

%some small error
total_input_error = sum(u1_tls-u2_tls)


% now estimate ML parameters (Tims trick)
phi1_est = inv(D1' * D1) * D1' * u1
fprintf('\n vs\n')
phi1
fprintf('\n and \n')
phi2_est = inv(D2' * D2) * D2' * u2
fprintf('\n vs\n')
phi2


%% Now find maximum discirminatability parameters phi1-phi2

%number of particles
n_pf = 1000;

% create a bunch of nearby phi particles
for ii = 1:length(phi1_est)
    phi1_pf(:,ii) = phi1_est(ii) + (rand(n_pf,1)-.5)*0.3*phi1_est(ii);
end
for ii = 1:length(phi2_est)
    phi2_pf(:,ii) = phi2_est(ii) + (rand(n_pf,1)-.5)*0.3*phi2_est(ii);
end



% %plot em
% figure(4)
% subplot(1,3,1);
% scatter(1:length(phi1_pf(:,1)),phi1_pf(:,1))
% hold on
% plot(1:length(phi1_pf(:,1)),ones(length(phi1_pf(:,1)),1)*phi1_est(1),'r')
% hold off
% title('phi 1, a0 particles')
% xlabel('pf #')
% ylabel('value')
% 
% subplot(1,3,2);
% scatter(1:length(phi1_pf(:,2)),phi1_pf(:,2))
% hold on
% plot(1:length(phi1_pf(:,2)),ones(length(phi1_pf(:,2)),1)*phi1_est(2),'r')
% hold off
% title('phi 1, a1 particles')
% xlabel('pf #')
% ylabel('value')
% 
% subplot(1,3,3);
% scatter(1:length(phi1_pf(:,3)),phi1_pf(:,3))
% hold on
% plot(1:length(phi1_pf(:,3)),ones(length(phi1_pf(:,3)),1)*phi1_est(3),'r')
% hold off
% title('phi 1, a2 particles')
% xlabel('pf #')
% ylabel('value')
% 
% figure(5)
% subplot(1,3,1);
% kk = 1;
% scatter(1:length(phi2_pf(:,kk)),phi2_pf(:,kk))
% hold on
% plot(1:length(phi2_pf(:,kk)),ones(length(phi2_pf(:,kk)),1)*phi2_est(kk),'r')
% hold off
% title('phi 2, a0 particles')
% xlabel('pf #')
% ylabel('value')
% 
% subplot(1,3,2);
% kk = 2;
% scatter(1:length(phi2_pf(:,kk)),phi2_pf(:,kk))
% hold on
% plot(1:length(phi2_pf(:,kk)),ones(length(phi2_pf(:,kk)),1)*phi2_est(kk),'r')
% hold off
% title('phi 2, a1 particles')
% xlabel('pf #')
% ylabel('value')
% 
% subplot(1,3,3);
% kk = 3;
% scatter(1:length(phi2_pf(:,kk)),phi2_pf(:,kk))
% hold on
% plot(1:length(phi2_pf(:,kk)),ones(length(phi2_pf(:,kk)),1)*phi2_est(kk),'r')
% hold off
% title('phi 2, a2 particles')
% xlabel('pf #')
% ylabel('value')


% now find which particle pair has the best discriminability
%particle filter way.
% max_delta = 0;
% max_index = 1;
% for jj = 1 :n_pf
%     max_test(jj) = sum( abs(0*D1*phi1_pf(jj,:)' - D1*phi2_pf(jj,:)') + abs(D2*phi2_pf(jj,:)'*0 - D2*phi1_pf(jj,:)') );
%     if(max_test(jj) >= max_delta)
%        max_delta = max_test(jj);
%        max_index = jj;
%     end
% end

%exhaustive way
max_delta = 0;
max_index = 1;
lambda = 1; %weighting
pcntd = 0.3; %limits
steps = 0.1; %steps
jj = 1;
for phi1_1 = (1-pcntd)*phi1_est(1):steps*phi1_est(1):(1+pcntd)*phi1_est(1)
    for phi1_2 = (1-pcntd)*phi1_est(2):steps*phi1_est(2):(1+pcntd)*phi1_est(2)
        for phi1_3 = (1-pcntd)*phi1_est(3):steps*phi1_est(3):(1+pcntd)*phi1_est(3)
            for phi2_1 = (1-pcntd)*phi2_est(1):steps*phi2_est(1):(1+pcntd)*phi2_est(1)
                for phi2_2 = (1-pcntd)*phi2_est(2):steps*phi2_est(2):(1+pcntd)*phi2_est(2)
                    for phi2_3 = (1-pcntd)*phi2_est(3):steps*phi2_est(3):(1+pcntd)*phi2_est(3)
                        phi1_ex(jj,:) = [phi1_1,phi1_2,phi1_3];
                        phi2_ex(jj,:) = [phi2_1,phi2_2,phi2_3];
                        max_test(jj) = sum( abs(D1*phi1_ex(jj,:)' - lambda*D1*phi2_ex(jj,:)') + abs(D2*phi2_ex(jj,:)' - lambda*D2*phi1_ex(jj,:)') );
                        if(max_test(jj) >= max_delta)
                           max_delta = max_test(jj);
                           max_index = jj;
                        end
                        jj = jj +1;
                    end
                end
            end
        end
    end
end

% output new parameters with maximum discriminative powers
phi1_md = phi1_ex(max_index,:)
fprintf('\n vs\n')
phi1
fprintf('\n and \n')
phi2_md = phi2_ex(max_index,:)
fprintf('\n vs\n')
phi2
max_delta

%plot em, exhaustive
figure(4)
subplot(1,3,1);
scatter(1:length(phi1_ex(:,1)),phi1_ex(:,1))
hold on
plot(1:length(phi1_ex(:,1)),ones(length(phi1_ex(:,1)),1)*phi1_est(1),'r')
hold on
plot(1:length(phi1_ex(:,1)),ones(length(phi1_ex(:,1)),1)*phi1_md(1),'g')
hold off
title('phi 1, a0 particles')
xlabel('pf #')
ylabel('value')

subplot(1,3,2);
scatter(1:length(phi1_ex(:,2)),phi1_ex(:,2))
hold on
plot(1:length(phi1_ex(:,2)),ones(length(phi1_ex(:,2)),1)*phi1_est(2),'r')
hold on
plot(1:length(phi1_ex(:,2)),ones(length(phi1_ex(:,2)),1)*phi1_md(2),'g')
hold off
title('phi 1, a1 particles')
xlabel('pf #')
ylabel('value')

subplot(1,3,3);
scatter(1:length(phi1_ex(:,3)),phi1_ex(:,3))
hold on
plot(1:length(phi1_ex(:,3)),ones(length(phi1_ex(:,3)),1)*phi1_est(3),'r')
hold on
plot(1:length(phi1_ex(:,3)),ones(length(phi1_ex(:,3)),1)*phi1_md(3),'g')
hold off
title('phi 1, a2 particles')
xlabel('pf #')
ylabel('value')

figure(5)
subplot(1,3,1);
kk = 1;
scatter(1:length(phi2_ex(:,kk)),phi2_ex(:,kk))
hold on
plot(1:length(phi2_ex(:,kk)),ones(length(phi2_ex(:,kk)),1)*phi2_est(kk),'r')
hold on
plot(1:length(phi2_ex(:,kk)),ones(length(phi2_ex(:,kk)),1)*phi2_md(kk),'g')
hold off
title('phi 2, a0 particles')
xlabel('pf #')
ylabel('value')

subplot(1,3,2);
kk = 2;
scatter(1:length(phi2_ex(:,kk)),phi2_ex(:,kk))
hold on
plot(1:length(phi2_ex(:,kk)),ones(length(phi2_ex(:,kk)),1)*phi2_est(kk),'r')
hold on
plot(1:length(phi2_ex(:,kk)),ones(length(phi2_ex(:,kk)),1)*phi2_md(kk),'g')
hold off
title('phi 2, a1 particles')
xlabel('pf #')
ylabel('value')

subplot(1,3,3);
kk = 3;
scatter(1:length(phi2_ex(:,kk)),phi2_ex(:,kk))
hold on
plot(1:length(phi2_ex(:,kk)),ones(length(phi2_ex(:,kk)),1)*phi2_est(kk),'r')
hold on
plot(1:length(phi2_ex(:,kk)),ones(length(phi2_ex(:,kk)),1)*phi2_md(kk),'g')
hold off
title('phi 2, a2 particles')
xlabel('pf #')
ylabel('value')


%assign difference vector
phi_d = phi1_md - phi2_md;

%% Run on some new data

%test by comparing with parameters close to either original
% phi1 = [100;30;1];
% phi2 = [15;200;1];
phi_new = [35;150;1];
phi = phi_new;

% SIMULINK
fprintf('Running simulation on new data...')

sim('SimpleModel1.slx');
fprintf(' DONE\n')

% no noise, get output stuff
u_new = input_out.Data(:,1); 
x_new = state.Data(:,3); 
xdot_new = state.Data(:,2); 
xdotdot_new = state.Data(:,1); 
state_new = state.Data;

%add in some noise
x_new = x_new + (rand(size(x_new))-.5)*0.01*max(x_new); 
xdot_new = xdot_new + (rand(size(xdot_new))-.5)*0.01*max(xdot_new);
xdotdot_new = xdotdot_new + (rand(size(xdotdot_new))-.5)*0.01*max(xdotdot_new);

Dx = [x_new,xdot_new,xdotdot_new];

% find the maximized delta
delta1 = sum((D1 - Dx)*phi_d');
delta2 = sum((Dx - D2)*phi_d');

figure(23)
plot(t, (D1 - Dx)*phi_d', 'r.'); hold on
plot(t, (Dx - D2)*phi_d', 'b.'); hold on
xlabel('t [s]')
ylabel('|Dx - Di| * Phi')
legend('D1-Dx', 'Dx-D2')



%compare deltas
if(delta1 > delta2)
   fprintf('Estimate is Class 1')
elseif(delta2 >= delta1)
    fprintf('Estimate is Class 2')
end