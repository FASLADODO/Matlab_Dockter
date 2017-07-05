% Homework 7
% 6.23
% Rod Dockter
% ME 5231

%% Part A)
K=1;
num=K;
den=[1 3 3 1];
sys=tf(num,den);
nyquist(sys);
%% Part B)
K=1;
num=K;
den=[1 1 0];
sys=tf(num,den);
nyquist(sys);
%% Part C)
K=1;
num=K;
den=[1 1 -1 -1];
sys=tf(num,den);
nyquist(sys);
%% Part D)
K=1;
num=K;
den=[1 1 -2];
sys=tf(num,den);
nyquist(sys);
%% Part E)
K=1;
num = K*[2 -6 -6 2];
den = [1 -3 3 -1];
sys=tf(num,den);
nyquist(sys);
%% Part F)
K=1;
num = conv([1 1],[1 0 3]);
den = K*[4 -12 12 -4];
sys=tf(num,den);
nyquist(sys);
%% Part G)
K=1;
num = K*[1 0 1];
den = [1 -3 3 -1];
sys=tf(num,den);
nyquist(sys);
%%
