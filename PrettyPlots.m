function [] = PrettyPlots(filename)

%eg:
%filename = 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\mrdPublications\ICRA_Letters_Grasper\Images\alpha_combined'

fontsize = 12;
FIG_W = 3.5;     % Width of actual figure  
FIG_H = 2;     % Height of actual figure
FIG_UNITS = 'inches'; % units for W&H
FIG_RES = 400; % figure resolution in dpi

% set it's W and H w/o messing up the position on the screen
set(gcf,'PaperPositionMode','auto', 'units', FIG_UNITS)
FIG_SZ = get(gcf, 'position');
FIG_SZ(3:end) = [FIG_W FIG_H];
set(gcf, 'position', FIG_SZ);

% Save the figure to file
if(nargin == 1)
    print(gcf, filename, '-dpng', ['-r' num2str(FIG_RES)]) % Raster
end

end