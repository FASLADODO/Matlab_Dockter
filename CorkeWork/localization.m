%% localization and mapping

p0 = diag([0.01, 0.01, 0.005].^2);
W = diag( [0.1, 1*pi/180].^2);
V = diag([0.02, 0.5*pi/180].^2);
map = Map(20);
veh = Vehicle(W);

veh.add_driver(RandomPath(map.dim));
sensor = RangeBearingSensor(veh,map,W);
ekf = EKF(veh,V,p0, sensor, W, []);

ekf.run(1000);

map.plot();
veh.plot_xy();
ekf.plot_map(5,'g');
ekf.plot_xy('r');

%% particle filter

p0 = diag([0.01, 0.01, 0.005].^2);
W = diag( [0.1, 1*pi/180].^2);
V = diag([0.005, 0.5*pi/180].^2);
map = Map(20);
veh = Vehicle(W);

veh.add_driver(RandomPath(10));
sensor = RangeBearingSensor(veh,map,V);

Q = diag([0.1, 0.1, 1*pi/180].^2);
L = diag([0.1, 0.1]);

pf = ParticleFilter(veh, sensor, Q, L, 1000);

pf.run(1000);

map.plot();
veh.plot_xy('b');
pf.plot_xy('r');

plot(pf.std(1:100,:));