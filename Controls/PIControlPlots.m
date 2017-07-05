%% overdamped
figure(1)

s= tf('s');

zero = 8;

D = ((s+zero)/s)

G = (1/((s+2)*(s+4)))

rlocus(G*D)

%% underdamped
figure(2)

s= tf('s');

zero = 8;

D = ((s+zero)/s)

G = (1/((s^2+2*s+6)))

rlocus(G*D)

%% critically damped
figure(3)

s= tf('s');

zero = 10;

D = ((s+zero)/s)

G = (1/((s+4)^2))

rlocus(G*D)

%% undamped
figure(4)

s= tf('s');

zero = 6;

D = ((s+zero)/s)

G = (1/((s^2+4)))

rlocus(G*D)