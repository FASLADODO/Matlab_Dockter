function K = Seongho(NGP,XCoord,YCoord,E,v,t)
x = XCoord;
y = YCoord;
D=(E/(1-v^2))*[1,v,0;v,1,0;0,0,(1-v)/2];
K=zeros(8);
B=zeros(3,8);

if NGP ==1
    GuassPoints=[0];
    Weights=[2];
end
if NGP ==2
    GuassPoints=[-0.57735027,0.57735027];
    Weights=[1,1];
end
if NGP ==3
    GuassPoints=[-0.77459666924148,0,0.77459666924148];
    Weights=[5/9,8/9,5/9];
end
if NGP ==4
    GuassPoints=[-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];
    Weights=[0.34785481,0.6521451549,0.6521451549,0.34785481];
end
i=1;
while  i<=NGP
    j=1;
    while j<=NGP
        s=GuassPoints(j);
        t=GuassPoints(i);
        W1=Weights(j);
        W2=Weights(i);
        a=0.25*(y(1)*(s-1)+y(2)*(-1-s)+y(3)*(1+s)+y(4)*(1-s));
        b=0.25*(y(1)*(t-1)+y(2)*(1-t)+y(3)*(1+t)+y(4)*(-1-t));
        c=0.25*(x(1)*(t-1)+x(2)*(1-t)+x(3)*(1+t)+x(4)*(-1-t));
        d=0.25*(x(1)*(s-1)+x(2)*(-1-s)+x(3)*(1+s)+x(4)*(1-s));
        J=0.125*x*[0,1-t,t-s,s-1;t-1,0,s+1,-s-t;s-t,-s-1,0,t+1;1-s,s+t,-t-1,0]*transpose(y);
        N1S=0.25*(t-1);
        N1T=0.25*(s-1);
        N2S=0.25*(1-t);
        N2T=0.25*(-1-s);
        N3S=0.25*(1+t);
        N3T=0.25*(1+s);
        N4S=0.25*(-1-t);
        N4T=0.25*(1-s);
        B1=[a*N1S-b*N1T,0;0,c*N1T-d*N1S;c*N1T-d*N1S,a*N1S-b*N1T];
        B2=[a*N2S-b*N2T,0;0,c*N2T-d*N2S;c*N2T-d*N2S,a*N2S-b*N2T];
        B3=[a*N3S-b*N3T,0;0,c*N3T-d*N3S;c*N3T-d*N3S,a*N3S-b*N3T];
        B4=[a*N4S-b*N4T,0;0,c*N4T-d*N4S;c*N4T-d*N4S,a*N4S-b*N4T];
        B=(1/J)*[B1 B2 B3 B4];
        K=K+(t*W1*W2*transpose(B)*D*B*J);
        j=j+1;
    end
   i=1+i;
end
K;