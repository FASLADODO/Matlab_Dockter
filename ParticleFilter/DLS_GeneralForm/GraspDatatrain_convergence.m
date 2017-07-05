clear all
close all

D_E_soft = load('SoftTissOnline.csv');

dt = 0.001;

Strain = D_E_soft(:,1);
Straindot = D_E_soft(:,2);
Straindotdot = D_E_soft(:,3);

strainDot_calc = smooth(diff(Strain)/dt ,7);
strainDot_calc = [strainDot_calc(1); strainDot_calc];

strainDotDot_calc = smooth(diff(strainDot_calc)/dt ,7);
strainDotDot_calc = [strainDotDot_calc(1); strainDotDot_calc];

figure
quiver(Strain,Straindot,Straindot,Straindotdot);

figure
subplot(3,1,1)
plot(1:length(Strain),Strain,'r.');
title('strain')
subplot(3,1,2)
plot(1:length(Straindot),Straindot,'r.');
hold on 
plot(1:length(strainDot_calc),strainDot_calc,'b.');
hold off
title('Straindot')
subplot(3,1,3)
plot(1:length(Straindotdot),Straindotdot,'r.');
hold on 
plot(1:length(strainDotDot_calc),strainDotDot_calc,'b.');
hold off
title('Straindotdot')

%%

for ii = 1:length(D_E_soft)
    if(D_E_soft(:,4) == 1)
        time(ii) =  D_E_soft(ii,1);
        Strain(ii) =  D_E_soft(ii,3);
        Straindot(ii) =  D_E_soft(ii,2);
    end
    
end

dt = 0.001;

strainDot_calc = smooth(diff(Strain)/dt ,7);
strainDot_calc = [strainDot_calc(1); strainDot_calc];

figure
plot(time,Strain,'r.');
title('strain')

figure
h1 = plot(time,Straindot,'r.');
hold on
h2 = plot(time,strainDot_calc,'b.');
hold off
title('straindots')

legend([h1,h2],'online','calced')


