function K =stiffnessSeongho(NGP)
K=zeros(18);
syms S T
X=[0,1.5,1.5,0,0.75,1.5,0.75,0,0.75]; % coordinate of X
Y=[0,0,4,4,0,2,4,2,2]; % coordinate of Y
GaussPoints=[-0.7746,-0.7746;0.7746,-0.7746;0.7746,0.7746;-0.7746,0.7746;0,-0.7746;0.7746,0;0,0.7746;-0.7746,0;0,0];
W=[5/9 5/9;5/9 5/9;5/9 5/9;5/9 5/9;8/9 5/9;5/9 8/9;8/9 5/9;5/9 8/9;8/9 8/9];
E = 10^7;
v = 0.3;
t = 0.2;
D=(E/(1-v^2))*[1,v,0;v,1,0;0,0,(1-v)/2];
N1=(1-S)*(1-T)/4; % shape function in terms of S,T
N2=(1+S)*(1-T)/4;
N3=(1+S)*(1+T)/4;
N4=(1-S)*(1+T)/4;
N5=(1-S^2)*(1-T)/2;
N6=(1+S)*(1-T^2)/2;
N7=(1-S^2)*(1+T)/2;
N8=(1-S)*(1-T^2)/2;
N9=(1-S^2)*(1-T^2);
N=[N1 N2 N3 N4 N5 N6 N7 N8 N9];
N0=[N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0 N9 0;0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0 N9];
X1=N1*X(1)+N2*X(2)+N3*X(3)+N4*X(4)+N5*X(5)+N6*X(6)+N7*X(7)+N8*X(8)+N9*X(9);
Y1=N1*Y(1)+N2*Y(2)+N3*Y(3)+N4*Y(4)+N5*Y(5)+N6*Y(6)+N7*Y(7)+N8*Y(8)+N9*Y(9);
J=[diff(X1,S), diff(Y1,S);diff(X1,T), diff(Y1,T)];
i=1;
while i <= NGP  
    D1=(1/det(J))*[diff(Y1,T)*diff(N(i),S)-diff(Y1,S)*diff(N(i),T),0;0,diff(X1,S)*diff(N(i),T)-diff(X1,T)*diff(N(i),S);diff(X1,S)*diff(N(i),T)-diff(X1,T)*diff(N(i),S),diff(Y1,T)*diff(N(i),S)-diff(Y1,S)*diff(N(i),T)];
    B=D1*N0;
    tempK=transpose(B)*D*B*t*W(i,1)*W(i,2)*det(J); 
    Stemp=GaussPoints(i,1);
    Ttemp=GaussPoints(i,2); 
    temp2K= subs(tempK,{S,T},{Stemp,Ttemp}); % call function stiffness
    K=K+temp2K;
    i=i+1;
end
K;