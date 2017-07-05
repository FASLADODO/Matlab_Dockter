function [K,alpha,T,mag] = lagcomp(w,factor)
K=w*sqrt((w/5)^2+1)*sqrt((w/50)^2+1);
alpha=(1/(0.01*K));
wcorn=w/factor;
T=1/wcorn;
num=alpha*sqrt((T*wcorn)^2+1);
den=sqrt((alpha*T*wcorn)^2+1);
mag=num/den;