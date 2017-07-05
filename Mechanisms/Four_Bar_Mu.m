%Four bar using erdmans mu*psi method

A0 = -495.000000 + -167.000000 *i ;
A1 = -406.915679 + -208.684691 *i ;
B0 = -257.000000 + -257.000000 *i ;
B1 = -343.761974 + -84.904776 *i ;
WA = 97.449787 ;
ZA = 440.659051 ;
WB = 124.769066 ;
ZB = 304.862868 ;
PP1 = -279 + 213*i;
PP2 = -88 + 163*i;
PP3 = -63 - 95*i;

%Specify this stuff...
L2 = 3;
L3 = 5;
L4 = 7;
L1 = 6;
theta1 = 0;
theta2 = 80;
elbow = 0; %1 for up, 0 for down
mu = 1;
prev_mu = 1;
index = 0;
dir = -1;

%Initiliaze
x1 = 0;
y1 = 0;
x3 = 0;
y3 = 0;
x4 = L1*cosd(theta1);
y4 = L1*sind(theta1);
x2 = L2*cosd(theta2);
y2 = L2*sind(theta2);


hf = figure('color','white');
while 1
    %Update Input
    theta2 = theta2 + dir*1;
    x2 = L2*cosd(theta2);
    y2 = L2*sind(theta2);
    
    %Update vector R7 (B0 to A)
    r7angle = atan2( ( y2 - y4) , (x2 - x4)) * (180/pi);
    r7mag = sqrt( (x2 - x4)^2 + ( y2 - y4)^2 );

    %Determine psi and theta4 angles
    psi = acosd((L4^2 + r7mag^2 - L3^2) / (2*L4*r7mag));
    if(psi > 180)
       psi = psi - 360;
    elseif(psi < -180)
       psi = psi + 360;
    end
    mu = sign(psi);
    
    if(index == 0)
       prev_mu = mu; 
    end
    
    %check for change in mu (singularity or circuit defect)
    if(prev_mu ~= mu)
        dir = dir * (-1); 
        theta2 = theta2 + dir*1;
        x2 = L2*cosd(theta2);
        y2 = L2*sind(theta2);

        %Update vector R7 (B0 to A)
        r7angle = atan2( ( y2 - y4) , (x2 - x4)) * (180/pi);
        r7mag = sqrt( (x2 - x4)^2 + ( y2 - y4)^2 );

        %Determine psi and theta4 angles
        psi = acosd((L4^2 + r7mag^2 - L3^2) / (2*L4*r7mag));
        if(psi > 180)
           psi = psi - 360;
        elseif(psi < -180)
           psi = psi + 360;
        end
        mu = sign(psi);
    end
    prev_mu = mu;
    theta4 = r7angle + psi;
  
    %Determine point B coordinates and theta3
    x3 = L4*cosd(theta4) + x4; %offset by ground 2
    y3 = L4*sind(theta4) + y4; %offset by ground 2
    theta3 = atand( (y3 - y2) / (x3 - x2));
    
    %Plot
    line([x1,x2],[y1,y2],'LineWidth',4,'Color',[0 0 0])
    hold on
    line([x2,x3],[y2,y3],'LineWidth',4,'Color',[0 0 1])
    hold on
    line([x4,x3],[y4,y3],'LineWidth',4,'Color',[0 1 0])
    hold off
    axis([-15, 15, -15, 15])
    title(['Branch is mu= ',num2str(mu)])

    %Animate
    refreshdata(hf,'caller')
    drawnow
    pause(0.001);
    clf

    index = index + 1;
end