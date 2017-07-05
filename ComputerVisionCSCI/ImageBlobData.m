function [meanim,minim,maxim,rangeim] = ImageBlobData(image)
%image is a binary image 0-255

    %image size
    [RR,CC] = size(image);

    BlobLocations = [];
    %loop through all image points
    for yy = 1:RR
        for xx = 1:CC
            temp = image(xx,yy);
            if(temp >= 255)
                %add to our array
                BlobLocations = [BlobLocations; xx,yy ];
            end
        end
    end
    
    %save image blob data
    minim = min(BlobLocations);
    maxim = max(BlobLocations);
    meanim = mean(BlobLocations);
    rangeim = range(BlobLocations);
end