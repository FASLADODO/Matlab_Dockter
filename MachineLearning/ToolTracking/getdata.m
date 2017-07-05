%Get Data in mat

D = load('data.txt');

T = D(:,1);
X = D(:,[2,3,4]);

figure
scatter(1:length(T),T);
ylabel('ms')
title('timestamps')

dt = mean(diff(T))
hz = 1000/dt

figure
scatter3(X(:,1),X(:,2),X(:,3));
xlabel('x(inches)')
ylabel('y(inches)')
zlabel('z(inches)')
title('raw tracking data')

span = 21;
XF1 = smooth(X(:,1),span);
XF2 = smooth(X(:,2),span);
XF3 = smooth(X(:,3),span);
XF = [XF1, XF2, XF3];

figure
scatter3(XF(:,1),XF(:,2),XF(:,3));
xlabel('x(inches)')
ylabel('y(inches)')
zlabel('z(inches)')
title('filter tracking data')

%%

%get spectrum 
[pxx,w] = periodogram(X,[],[],hz);
figure
plot(w,10*log10(pxx))
% plot((w/pi)*hz,pxx)
xlabel('Hz')
ylabel('dB')
title('Periodogram of tool motion data')
legend('X','Y','Z')

%% Filters

% Butterworth - maximize flatness of pass band
% Chebychev - maximize sharpness of cutoff




