% test kl divergence


p1 = 0:0.01:1;

p2 = 0:0.01:1;

[XX,YY] = meshgrid(p1,p2);

D = [XX(:),YY(:)];

pwithin = D(:,1);
pbetween = D(:,2);

Diff1 = log10(pwithin./(pbetween)); %LOG LIKELIHOOD
Diff2 = pwithin.*log10(pwithin./(pbetween)); %REGULAR KL
Diff3 = (pwithin+pbetween).*log10(pwithin./(pbetween)); %NOT ACTUALLY KL 
Diff4 = pwithin.*log10(pwithin./pbetween) + pbetween.*log10(pbetween./pwithin); %SYMMETRIC KL

figure
Surface3D(D(:,1),D(:,2),Diff1);
xlabel('x1')
ylabel('x2')
zlabel('metric')
title('KL metric 1')

figure
Surface3D(D(:,1),D(:,2),Diff2);
xlabel('x1')
ylabel('x2')
zlabel('metric')
title('KL metric 2')

figure
Surface3D(D(:,1),D(:,2),Diff3);
xlabel('x1')
ylabel('x2')
zlabel('metric')
title('KL metric 3')

figure
Surface3D(D(:,1),D(:,2),Diff4);
xlabel('x1')
ylabel('x2')
zlabel('metric')
title('KL metric 4')

ComputeRBFDifference(0.4,0.2)

ComputeRBFDifference(0.2,0.1)