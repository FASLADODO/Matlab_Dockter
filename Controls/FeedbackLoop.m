function [y, determ, invy] = FeedbackLoop(s)
syms s
G1=[(1-s)/s 1/(s+1);1/s -s/(s+1)];
G2=[1 0;0 1];
y=eye(2)+G1*G2;
determ=det(y);
invy=inv(eye(2)+G1*G2);


