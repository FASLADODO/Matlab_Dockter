function K =Triangle(XCoord,YCoord,Stress,Strain,v,t,A)
%%Present coordinate in counterclock wise fashion
%%enter 0 or 1 depending if you want stress of strain
x=XCoord;
y=YCoord;
D=zeros(3);
if Stress==1
    D=(E/(1-v^2))*[1,v,0;v,1,0;0,0,(1-v)/2];
end
if Strain==1
   D=(E/((1+v)*(1-2*v)))*[1-v,v,0;v,1-v,0;0,0,(1-2*v)/2]; 
end
Bi=y(2)-y(3);
Bj=y(3)-y(1);
Bm=y(1)-y(2);
Gi=x(3)-x(2);
Gj=x(1)-x(3);
Gm=x(2)-x(1);
B=(1/(2*A))*[Bi,0,Bj,0,Bm,0;0,Gi,0,Gj,0,Gm;Gi,Bi,Gj,Bj,Gm,Bm];
K=t*A*transpose(B)*D*B;