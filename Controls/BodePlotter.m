function [wbw,wbwnp] = BodePlotter()
%6.12 Plots bode daigram for Transfer Function with additional poles
%for 5 different values of p specified
p=[0.01,0.1,1,10,100];
[wbw]=zeros(5);
for i= 1:5
    sys(i)=tf(1,[1/p(i),1+(1/p(i)),1+(1/p(i)),1]);
    wbw(i)=bandwidth(sys(i));
end
bodeplot(sys(1),sys(2),sys(3),sys(4),sys(5));
sys2=tf([0,0,1],[1,1,1]);
wbwnp=bandwidth(sys2)
