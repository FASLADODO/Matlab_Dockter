function [zoomim] = ZoomImageCenter(image,scale)
%will scale the image by the x and y parameters in scale
%then filters the scaled image
%then crops back the scaled image to be the same size as the orginal
%just zoomed into the center

    %input image size
    [RR,CC,depth] =size(image);
    
    %new scale
    newcols = round(scale(1)*CC)
    newrows = round(scale(2)*RR)
    zoomim = zeros(newrows, newcols, depth, class(image));
    
    %loop through all pixels in new picture
    for cc = 1:newcols
       for rr = 1:newrows 
           %figure where this pixel is coming from in old picture
            oldcoordc = round(cc/scale(1));
            oldcoordr = round(rr/scale(2));
            %just in case we go outside
            if(oldcoordc < 1)
                oldcoordc = 1;
            elseif(oldcoordc > CC)
                oldcoordc = CC;
            end
            if(oldcoordr < 1)
                oldcoordr = 1;
            elseif(oldcoordr > RR)
                oldcoordr = RR;
            end
            
            %stash old pixels via interpolation
            zoomim(rr,cc,:) = image(oldcoordr,oldcoordc,:);
       end
    end
    
    %anti-aliasing filter
    for dd = 1:depth;
        tempim = zoomim(:,:,dd);
        tempsmooth = medfilt2(tempim);
        zoomim(:,:,dd) = tempsmooth;
    end
    
    %now crop back to original size
    centercol = round(newcols / 2)
    centerrow = round(newrows / 2)
    if(newcols >= CC)
        halfcol = round(CC /2);
    else
        halfcol = round(newcols /2);
    end
    if(newrows >= RR)
        halfrow = round(RR /2);
    else
        halfrow = round(newrows /2);
    end
    %figure out our new top left 
    newtopleft = [centerrow-halfrow, centercol-halfcol];
    if(newtopleft(1) < 1)
        newtopleft(1) = 1;
    end
    if(newtopleft(2) < 1)
        newtopleft(2) = 1;
    end
    %make it exatcly the same size
    newidxrow = [ newtopleft(1):(newtopleft(1) + RR - 1) ];
    newidxcol = [ newtopleft(2):(newtopleft(2) + CC - 1) ];

    %only keep the orginal image size in the center of scaled image
    zoomim = zoomim( newidxrow, newidxcol, :);
    

end






