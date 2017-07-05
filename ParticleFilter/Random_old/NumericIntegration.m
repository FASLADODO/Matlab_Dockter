function[K] = NumericIntegration (NGP,XCoord,YCoord,E,v,t,Weights,Gpoints)
%%Enter X coordinate and Y coordinate in increasing order array
%%Weights and Gauss Points are optional inputs
W1=0;
W2=0;
X=XCoord;
Y=YCoord;
%%Spoints for Xi points
%%Npoints for Eta points
if nargin == 6
    switch NGP
        case 1
            Npoints=[0];
            Spoints=[0];
            W=[2];
        case 2
            Npoints=[-0.57735027918962,0.57735027918962];
            Spoints=[-0.57735027918962,0.57735027918962];
            W=[1,1];
        case 3
            Npoints=[-0.77459666924148,0,0.77459666924148];
            Spoints=[-0.77459666924148,0,0.77459666924148];
            W=[5/9,8/9,5/9];
        case 4
            Npoints=[-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];
            Spoints=[-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];
            W=[0.34785481,0.6521451549,0.6521451549,0.34785481];
    end
end
if nargin == 8
   Npoints=Gpoints;
   Spoints=Gpoints;
   W=Weights;
end
D=(E/(1-v^2))*[1,v,0;v,1,0;0,0,(1-v)/2];
i=1;
B=zeros(3,8);
K=zeros(8);
while i <= NGP
    j=1;
    while j <= NGP
        s=Spoints(j);
        n=Npoints(i);
        W1=W(j);
        W2=W(i);
        j11=(1/4)*((n-1)*X(1)+(1-n)*X(2)+(n+1)*X(3)+(-n-1)*X(4));
        j12=(1/4)*((n-1)*Y(1)+(1-n)*Y(2)+(n+1)*Y(3)+(-n-1)*Y(4));
        j21=(1/4)*((s-1)*X(1)+(-s-1)*X(2)+(s+1)*X(3)+(1-s)*X(4));
        j22=(1/4)*((s-1)*Y(1)+(-s-1)*Y(2)+(s+1)*Y(3)+(1-s)*Y(4));
        J=[j11,j12;j21,j22];
        Jdet=det(J);
        B1=[j22*(n-1)-j12*(s-1);0;j11*(s-1)-j21*(n-1)];
        B2=[0;j11*(s-1)-j21*(n-1);j22*(n-1)-j12*(s-1)];
        B3=[j22*(1-n)-j12*(-s-1);0;j11*(-s-1)-j21*(1-n)];
        B4=[0;j11*(-s-1)-j21*(1-n);j22*(1-n)-j12*(-s-1)];
        B5=[j22*(n+1)-j12*(s+1);0;j11*(s+1)-j21*(n+1)];
        B6=[0;j11*(s+1)-j21*(n+1);j22*(n+1)-j12*(s+1)];
        B7=[j22*(-n-1)-j12*(1-s);0;j11*(1-s)-j21*(-n-1)];
        B8=[0;j11*(1-s)-j21*(-n-1);j22*(-n-1)-j12*(1-s)];
        B=(1/4)*(1/Jdet)*[B1,B2,B3,B4,B5,B6,B7,B8];
        Btran=transpose(B);
        K=K+(t*W1*W2*Btran*D*B*Jdet);
        j=j+1;
    end
    i=i+1;
end
K;
end

