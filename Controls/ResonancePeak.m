function[Mr]=ResonancePeak()
%6.11 Frequency Response
[a]=[0.01,0.1,1,10,100];
w=logspace(-1,1);
[T]=zeros(5);
[Mr]=zeros(5);
[J]=zeros(5);
for i = 1:5
    sys = tf([1/a(i),1],[1,1,1]);
    y=bode(sys,w);
    [J(i),T(i)]=max(y);
end
for k=1:5
    Mr(k)=J(k)-1.0;
end
