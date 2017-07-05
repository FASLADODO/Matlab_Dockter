% Sensitivity for closed loop transfer function

syms s a k

%closed loop transfer function
T = k/(s^2+a*s+k);

%dT/da
deriv = diff(T,a);

% sensivity per eq 7.75 in nise
sensitivity = (a/T)*deriv
