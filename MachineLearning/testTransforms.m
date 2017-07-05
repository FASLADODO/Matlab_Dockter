P = [20,20;
    20,40;
    40,40;
    40,20];



T = Transformation2D(0.1,[2,2],1)

% T = [1,0,0;
%     0,1,0;
%     0,0,1];

Pprime = Transform2x2(P,T);


figure
plot(P(:,1),P(:,2),'r*')
hold on
plot(Pprime(:,1),Pprime(:,2),'b*')
hold off
