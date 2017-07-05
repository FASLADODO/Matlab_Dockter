% function toTexFig(figName,captionText, width, rescale) 
%   This function will export .eps and .jpg files of the current figure named
%   figName;  It will also spit out text to stdout which is LaTeX code for
%   the fig...and captionText if it is provided, This includes:
%
% %%%%%% FIGURE:  figName %%%%%%%%%%%%%%
%   \begin{figure}[htbp]			 % look into how the position of figures are determined
%   \centering
%       \includegraphics[]{figName}  % Should take eps as the default
%       \caption{captionText}
%       \label(fig:figName)
%   \end{figure}

function toTexFig(figName,captionText, width, rescale) %, LegendText
Res = '-r600';  % default figure resolution
% Extract svn root dir
d=strrep(pwd,'\','/');
figDir = [d(1:strfind(d,'/surgenome/')) 'surgenome/dissertation/fig/'];

%lh=legend('PegTx','Cutting','Suturing');
%set(lh,'orientation', 'horizontal');
%set(lh,'position', [0.2296    0.9372    0.6095    0.0443]);
% set(lh, 'FontSize', 14);
%set(gcf,'PaperPositionMode','auto');
%print('-dpdf',Res, [figDir figName])
% print('-dpsc','-loose', [figDir figName]) %ps file, color
% '-painters'
% '-zbuffer'
% '-opengl'
print('-depsc','-loose', [figDir figName]) %eps file, color
%print('-depsc','-loose', [figDir figName '-loose']) 
% paper width should be standard constant: 
if (rescale)
    if length(width)==1
        w = width; % inches        
        figPos = get(gcf,'position');
        h = figPos(4)/figPos(3) *w; % aspect ratio of fig mapped to paper
    else 
        w=width(1);
        h=width(2);
    end
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'paperposition',[ 0 0  w h ])
    % get(gcf, 'papertype')
    set(gcf,'papersize', [h w]);
end
% get(gcf, 'papertype')
print('-dpng',Res, [figDir figName])
%print('-dpdf',Res, [figDir figName])
%%
fprintf(1,['\n\n\n'...
    '%%%%%%%%%%%%%%%%%%%%%%%% FIGURE:  %s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'...
    '\\begin{figure}[htbp]   \n'...
    '   \\centering \n'...
    '   \\includegraphics[width=%gin]{%s}  %% Should take eps as the default \n'...
    '   \\caption{%s} \n'...
    '   \\label{fig:%s} \n' ...
    '\\end{figure} \n'...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'...
    '\n\n'],figName, width, figName, captionText, figName);

