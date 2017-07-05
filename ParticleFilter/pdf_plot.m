function [ ] = pdf_plot(pdfMatrix,x,y)

hold on;grid on;

% Contour of probability distribution
[C,h] = contourf(x,y,pdfMatrix);
% colorbar;

% Contour of normal acceleration
acenMatrix = repmat(x.^2,length(y),1).*10.^(repmat(y',1,length(x)));
[C1,h1] = contour(x,y,acenMatrix,[1,2,4,8,16,32],'y-.','LineWidth',2);
clabel(C1,h1,'Color','y');

xlabel('Speed [cm/s]');
ylabel('Curvature [log_{10}(cm^{-1})]');
text(4,1.6,'Normal Acc [cm/s^2]','Color','r');