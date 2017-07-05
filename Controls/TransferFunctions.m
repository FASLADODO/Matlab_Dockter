%Transfer function stuff

s = tf('s');  % defines the complex Laplace domain variable s

% You can use it to form Transfer functions more intuitively:

G = 10/ (s^2 + 2*s + 1);  % My plant

% Then do one of these:

pzmap(G)

impulse(G)

step(G)

step(G/s) % ramp response

ltiview(G)  % utility that lets you get all you want (settling time, etc). 