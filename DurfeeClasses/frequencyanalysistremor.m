
% taking L1 norm of the 3 axis of acceleration
a = ax.*ax + ay.*ay + az.*az;

Time = (count - count(1)) * 0.01;


%%
 
plot(Time,a)
axis([2,9,0,3*10^5])
%found to be about 12 cycles per 200 time units (0.06 hz)
%%
fdatool
