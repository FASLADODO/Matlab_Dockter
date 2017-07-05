syms s
expand((s+1)*(s+3)*(s+50))

%%

% p s^2+6*s+10.2225
p = [1,3,7,5];
roots(p)

%%

k = 4

num = [k*2];
den = [1,2+k*2];

sys = tf(num,den);

w = 1.41;
[gain,phase ] = bode(sys,w)

bode(sys)

[Gm,Pm,Wgm,Wpm] = margin(sys)

wb = bandwidth(sys)


%%

num = [2];
den = [1,2];

sys = tf(num,den);

nyquist(sys)


%%
num = [3];
den = [1,2,-6];

sys = tf(num,den);

rlocus(sys)

%%
k=1000;
num = [k];
den = [1,3+k];

sys = tf(num,den);

step(sys)
%axis([-10,0,-10,10])

S = stepinfo(sys)

