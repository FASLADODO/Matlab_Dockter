% simulinkin driving

%sl_drivepoint

xg = [5,5]; %goal coordinate
x0 = [8,5,pi/2]; %initial pose

r = sim('sl_drivepoint')

path = r.find('yout');

plot( path(:,1), path(:,2))

%% simulink fly

mdl_quadrotor

sim('sl_quadrotor')

about(result)

plot(result(:,1), result(:,2))
