% Plotter for 4 bar mechanism


% setup
clear all
deg2rad = (pi/180);
rad2deg = (180/pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify Beta (angle between input positions) (deg)
beta2 = 36;
beta3 = 120; 

% Precision positions real+i*im (x,y) (mm)
PP1 = -279 + 213*i;
PP2 = -88 + 163*i;
PP3 = -63 - 95*i;

% ground pivots  (mm) complex vector form
A0 = -495.000000 + -167.000000 *i ;
B0 = -257.000000 + -257.000000 *i ;

%moving pivots
A1 = -406.915679 + -208.684691 *i;
B1 = -372.242811 + -76.220642 *i;

% link lengths
WA = 88.084320799021190 - 41.684690823813135*i;
ZA = 1.279156792009793e+02 + 4.216846908238136e+02*i;
WB = -1.152428107344502e+02 + 1.807793578765912e+02*i;
ZB = 9.324281073445037e+01 + 2.892206421234089e+02 * i;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%%%%

% get coupler length between moving pivots
C = abs(B1 - A1);

% ground length
G = abs(B0 - A0);

%initial input angle 
beta1 = angle(A1 - A0) * rad2deg;



hf = figure('color','white');
while 1
   for arg = beta1 : 1 : beta3
       [alpha, psi, Ax, Ay, Bx, By] = PositionAnalysis4Bar( G, abs(WA), C, abs(WB), arg );
       
        %new moving pivots
        A1 = A0 + WA *exp(i*arg*deg2rad);
        B1 = B0 + WB *exp(i*arg*deg2rad);
        
        %new coupler vectors
        CA = ZA*exp(i*alpha*deg2rad) + A1;
        CB = ZB*exp(i*alpha*deg2rad) + B1;
       
       %%%%% Plotting %%%%%%
        plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
        hold on
        plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
        hold on
        plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot(real(A0),imag(A0),'p','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
        hold on
        plot(real(B0),imag(B0),'p','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot(real(A1),imag(A1),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
        hold on
        plot(real(B1),imag(B1),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot([real(A0), real(A1)],[imag(A0), imag(A1)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(B0), real(B1)],[imag(B0), imag(B1)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(B1), real(CB)],[imag(B1), imag(CB)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(A1), real(CA)],[imag(A1), imag(CA)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(A1), real(B1)],[imag(A1), imag(B1)],'-','LineWidth',2,'Color',[0 0 0]);
        hold off
        

        axis([-800, 800, -800, 800])
        %Animate
        refreshdata(hf,'caller')
        drawnow
        pause(0.001);
        clf
   end
    
   for arg = beta3 : -1 : beta1
       [alpha, psi, Ax, Ay, Bx, By] = PositionAnalysis4Bar( G, WA, C, WB, arg );
       
        %new moving pivots
        A1 = A0 + WA *exp(i*arg*deg2rad);
        B1 = B0 + WB *exp(i*arg*deg2rad);
        
        %new coupler vectors
        CA = ZA*exp(i*alpha*deg2rad) + A1;
        CB = ZB*exp(i*alpha*deg2rad) + B1;
       
       %%%%% Plotting %%%%%%
        plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
        hold on
        plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
        hold on
        plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot(real(A0),imag(A0),'p','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
        hold on
        plot(real(B0),imag(B0),'p','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot(real(A1),imag(A1),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
        hold on
        plot(real(B1),imag(B1),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot([real(A0), real(A1)],[imag(A0), imag(A1)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(B0), real(B1)],[imag(B0), imag(B1)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(B1), real(CB)],[imag(B1), imag(CB)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(A1), real(CA)],[imag(A1), imag(CA)],'-','LineWidth',2,'Color',[0 0 0]);
        hold on
        plot([real(A1), real(B1)],[imag(A1), imag(B1)],'-','LineWidth',2,'Color',[0 0 0]);
        hold off
       
        axis([-800, 800, -800, 800])
        %Animate
        refreshdata(hf,'caller')
        drawnow
        pause(0.001);
        clf
       
   end
    
end




