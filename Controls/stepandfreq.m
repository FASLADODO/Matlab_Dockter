function [Mp,T] = stepandfreq()
%6.11 Step response
[a]=[0.01,0.1,1,10,100];
[T]=zeros(5);
[Mp]=zeros(5);
[J]=zeros(5);
for i = 1:5
    sys = tf([1/a(i),1],[1,1,1]);
    y=step(sys);
    [J(i),T(i)]=max(y);
end
for k=1:5
    Mp(k)=J(k)-1.0;
end

    