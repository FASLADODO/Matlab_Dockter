% Make up some landmarks then localize
close all
clear all

nlm = 50; %number of landmarks

%set initial landmarks
landmarkloc = [randn(nlm,1)*10, randn(nlm,1)*10 ];

thresh = min(range(landmarkloc))*0.1;
windowsz = min(range(landmarkloc))*0.4;

%set initial robot locations
%thetan = (2*pi).*rand(1,1);
%robloc = [randn(1,1)*10, randn(1,1)*10 , thetan];
robloc = [-20 + rand*2, -20 + rand*2, 0.9 + rand*0.1];


lochist_actual = [];
lochist_odom = [];
lochist_corrected = [];

figure
%landmarks
h1 = scatter(landmarkloc(:,1), landmarkloc(:,2), 30, 'r*');
title('OG Landmarks')

%%

ff = figure(2);

while(1)

    %Odom robot position change
    roblocNew_Actual = robloc + [rand*3,rand*3,randn*0.1 ];
    roblocNew = roblocNew_Actual + [rand*0.2,rand*0.2,randn*0.01]; %noise added

    windowlandmarks = DataInWindow(landmarkloc,roblocNew([1,2]),windowsz);

    [nlw,~] = size(windowlandmarks);

    %Now get noisy measurement signals (this would normally just come form the
    %laser scan
    for ii = 1:nlw
       thtemp =  atan2(windowlandmarks(ii,2) - roblocNew_Actual(2), windowlandmarks(ii,1) - roblocNew_Actual(1) );
       bear = thtemp - roblocNew_Actual(3) + rand*0.01;
       dist = sqrt( (windowlandmarks(ii,2) - roblocNew_Actual(2))^2 + (windowlandmarks(ii,1) - roblocNew_Actual(1))^2 ) + rand*0.1;
       measuremarks(ii,:) = [dist,bear];
    end

    [Update, CurrentLandmarks] = localization2D(roblocNew, landmarkloc, measuremarks,thresh);

    lochist_actual = [ lochist_actual; roblocNew_Actual];
    lochist_odom = [ lochist_odom; roblocNew];
    lochist_corrected = [ lochist_corrected; Update];

    disp('We good?')
    Accuracy = (Update - roblocNew_Actual) ./ roblocNew_Actual

    clf(ff)

    
    %landmarks
    h1 = scatter(landmarkloc(:,1), landmarkloc(:,2), 30, 'r*');
    hold on
    h2 = scatter(CurrentLandmarks(:,1), CurrentLandmarks(:,2), 30, 'b*');
    hold on
    %Start Points
    h3 = quiver(robloc(1),robloc(2),cos(robloc(3)),sin(robloc(3)),'k','LineWidth',3);
    hold on

    %actual location
    h4 = quiver(roblocNew_Actual(1),roblocNew_Actual(2),cos(roblocNew_Actual(3)),sin(roblocNew_Actual(3)),'g','LineWidth',3);
    hold on
    h4 = plot(lochist_actual(:,1), lochist_actual(:,2), 'g-.');
    hold on
    %odom location
    h5 = quiver(roblocNew(1),roblocNew(2),cos(roblocNew(3)),sin(roblocNew(3)),'m','LineWidth',3);
    hold on
    h5 = plot(lochist_odom(:,1), lochist_odom(:,2), 'm-.');
    hold on
    %corrected location
    h6 = quiver(Update(1),Update(2),cos(Update(3)),sin(Update(3)),'c','LineWidth',3);
    hold on
    h6 = plot(lochist_corrected(:,1), lochist_corrected(:,2), 'c-.');
    hold off
    legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'True Landmark','Measured landmarks','Previous Location','Location Actual','Odometry Location','Corrected Location');

    robloc = Update;

    w = waitforbuttonpress;
    
    if w == 0
        %keep going
    else
        break;
    end
end