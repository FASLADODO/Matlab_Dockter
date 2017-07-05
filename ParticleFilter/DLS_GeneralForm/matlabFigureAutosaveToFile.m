%% Demonstrates how to export your matlab figures to a file programatically
%  Really nice for LaTeX papers & dissertations.  

% The width and height of the figure that will appear in your paper.  i.e.,
% it will measure exactly 3"x2" on the page and at that size the font will 
% scale correctly: 12pt font will match 12pt font in text.
FIG_W = 3.5;     % Width of actual figure  
FIG_H = 1.8;     % Height of actual figure
FIG_UNITS = 'inches'; % units for W&H
FIG_RES = 600; % figure resolution in dpi

figure(123)

% Make a plot
plot(rand(3,100)'); xlabel('index'); ylabel('noise')
title('\bfMeaningless Data')

% set it's W and H w/o messing up the position on the screen
set(gcf,'PaperPositionMode','auto', 'units', FIG_UNITS)
FIG_SZ = get(gcf, 'position');
FIG_SZ(3:end) = [FIG_W FIG_H];
set(gcf, 'position', FIG_SZ);


% Save the figure to file
print(gcf, 'MyFigFilename', '-dpng', ['r' num2str(FIG_RES)]) % Raster
print(gcf, 'MyFigFilename', '-dsvg', ['r' num2str(FIG_RES)]) % Vector