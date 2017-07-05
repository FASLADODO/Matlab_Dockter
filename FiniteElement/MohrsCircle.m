%Rods Mohr Circle function
%Oct 2014
function [sigma1,sigma2,sigma3,Tmax,thetaP] = MohrsCircle(sigmaX, sigmaY, Txy)
    %circle values
    center = (sigmaX + sigmaY)/2;
    R = sqrt(((sigmaX-sigmaY)/2)^2+Txy^2);
    thetaP = (1/2) * atand((2*Txy)/(sigmaX-sigmaY));
    
    %compute principals
    sigma3 = 0;
    sigma1 = center + R;
    sigma2 = center - R;
    
    if(sigma2 > 0) %3D mohrs
        %recompute principals
        sigma3 = sigma2; 
        Tmax = sigma1/2;

        %Plot everything
        figure(1)
        circle(center,0,R,[0 0 0]);
        hold on 
        circle(sigma1/2,0,sigma1/2,[1 0.5 1]);
        hold on 
        circle(sigma3/2,0,sigma3/2,[1 0 0.5]);
        hold on 
        plot(center,0,'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
        hold on
        plot(sigma1,0,'o','MarkerSize',5,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
        hold on
        plot(sigma2,0,'o','MarkerSize',5,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot(sigma1/2,Tmax,'o','MarkerSize',5,'MarkerFaceColor',[0 1 1],'Color',[0 1 1]);
        hold on
        plot([sigmaX,sigmaY],[Txy,-Txy],'-','Color',[0 0 0]);
        hold off

        %format
        grid on
        axis square
        legend('mohrs 1', 'mohrs total', 'mohrs 3' ,'center', 'sigma 1', 'sigma 2', 'T max')
        title('mohrs circle plot')
        xlabel('Normal Stress')
        ylabel('- Shear Stress')
        
    else %2D mohrs
        %shear max
        Tmax = R;

        %plot everything
        figure(1)
        circle(center,0,R,[1 0.5 0]);
        hold on 
        plot(center,0,'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
        hold on
        plot(sigma1,0,'o','MarkerSize',5,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
        hold on
        plot(sigma2,0,'o','MarkerSize',5,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
        hold on
        plot(center,Tmax,'o','MarkerSize',5,'MarkerFaceColor',[0 1 1],'Color',[0 1 1]);
        hold on
        plot([sigmaX,sigmaY],[Txy,-Txy],'-','Color',[0 0 0]);
        hold off

        %format
        grid on
        axis square
        legend('mohrs circle', 'center', 'sigma 1', 'sigma 2', 'T max')
        title('mohrs circle plot')
        xlabel('Normal Stress')
        ylabel('- Shear Stress')
    end

end

%helper function for plotting circles
function circle(x,y,r,color)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,'-','LineWidth',2,'Color', color);
end