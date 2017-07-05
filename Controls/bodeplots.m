syms s
expand((1/10)*s^2+(1/2)*s+1)


%%

num = [1,5];
den = [1,6,8];
sys = tf(num,den)

w = 10;

[mag,phase] = bode(sys,w)

bode(sys)

%%

num = [1,5];
den = [1,0,0];
sys = tf(num,den)
nyquist(sys)

%%

num = [400];
den = [1,10,0];
sys = tf(num,den)

rlocus(sys)

