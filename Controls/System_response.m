%Find System Specs, second order system

%T = wn^2 / (s^2 + 2*z*wn*s + wn^2)

%Enter denominator,  sys = [1, 2*z*wn, wn^2]

format long g

num = [1];
den = [1, 4, 10];

wn = sqrt(den(3));

z = den(2)/(2*wn);


%Display
natural_freq = wn
damping_ratio = z

rise_time_1 = 2.2 / (z*wn)
rise_time_2 = (1.76*z^3 - 0.417*z^2 + 1.039*z + 1) %fancy
%https://courses.engr.illinois.edu/ece486/lab/estimates/estimates.html

settle_time = 4 / (z*wn)

peak_time = pi / (wn*sqrt(1 - (z^2)))

percent_os = exp(-(z*pi / (sqrt(1 - (z^2))))) * 100

systf = tf(num,den)
step(systf)

%%
%enter percent overshoot & Tp, get out damping ratio
OS = 20;
T_p = 4;

z = (-log(OS/100))/(sqrt(pi^2+log(OS/100)^2))

w_n = pi/(T_p*sqrt(1-z^2))

%%

%Graph
%systf = tf([num],[den])

num = [0.04];
den = [1,0.02,0.04]
systf = tf(num,den)
step(systf)