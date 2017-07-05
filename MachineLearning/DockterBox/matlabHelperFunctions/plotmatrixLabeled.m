function [H,AX,BigAx,P,PAx] = plotmatrixLabeled(X,Title, Lbl)
[H,AX,BigAx,P,PAx]= plotmatrix(X);

% for each plot, determine a label
d=min(size(X));
for i=1:d.^2
    lbl{i}='';
end

top=d+[1:d] ;
side=[1:d]*d+1 ;
for i=1:length(top)
   lbl{top(i)}  = Lbl{i};
   lbl{side(i)} = Lbl{i};
end

% Find the handles to all axes. Refer to Axes Properties Section of the
% documentation for more detials on what each properties do.
axes_handles = findobj(gcf,'Type','Axes');
axes_posns = get(axes_handles(:),'Position');
axes_posns_mat = cell2mat(axes_posns);
xposns = sort(axes_posns_mat(:,1));
yposns = sort(axes_posns_mat(:,2));
xreq = [];
yreq = [];

% find the unique X-positions for the axes Position vector
for i =1:length(xposns)
    p = find(xposns(i)==xposns(:));
    if length(p)>1 && xposns(i) ~= NaN
        xreq = [xreq;xposns(i)];
        xposns(p) = NaN;
    end
end

% find the unique Y-positions
for i =1:length(yposns)
    p = find(yposns(i)==yposns(:));
    if length(p)>1 && yposns(i) ~= NaN
        yreq = [yreq;yposns(i)];
        yposns(p) = NaN;
    end
end

% Insert label or title
for i = 1:length(axes_handles)
%     ii=(i-1)/nScores;
%     jj=1+(ii-1)*nScores;
    p = get(axes_handles(i),'Position');
    if (p(1) == xreq(1)) && ~isempty(p(2) == yreq(:))
        set(gcf,'CurrentAxes',axes_handles(i)); ;
        ylabel(lbl(i))
        %if ii~=0;
            %ylabel(['#' num2str(ii)]);
            %ylabel(vidRevLabels{cols{c}(ii)+3})
            
        %end
        %title(vidRevLabels{i})
    end
    if (p(2) == yreq(end)) && ~isempty(p(1) == xreq(:))
        set(gcf,'CurrentAxes',axes_handles(i)); ;
        title(lbl(i))
        %if ii~=0;
            %title(['#' num2str(jj)]);
            %title(vidRevLabels{cols{c}(jj)+3})
            
        %end

    end
    if (p(1) == xreq(1)) && (p(2)== yreq(end))
        q = get(get(axes_handles(i),'children'),'marker');
        if q == '*'
            set(axes_handles(i),'visible','off');
        end
    end
end

suptitle(Title);