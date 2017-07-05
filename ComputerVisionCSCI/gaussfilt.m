% GAUSSFILT -  Small wrapper function for convenient Gaussian filtering
%
% Usage:  smim = gaussfilt(im, sigma)
%
% Arguments:   im - Image to be smoothed.
%           sigma - Standard deviation of Gaussian filter.
%
% Returns:   smim - Smoothed image.
%
% If called with sigma = 0 the function immediately returns with im assigned
% to smim


function smim = gaussfilt(im, sigma)
 
    if sigma < eps
        smim = im;
        return;
    end
    
    % If needed convert im to double
    if ~strcmp(class(im),'double')
        im = double(im);  
    end
    
    sze = max(ceil(6*sigma), 1);
    if ~mod(sze,2)    % Ensure filter size is odd
        sze = sze+1;
    end
    
    h = fspecial('gaussian', [sze sze], sigma);

    % Apply filter to all image channels
    smim = zeros(size(im));
    for n = 1:size(im,3)
        smim(:,:,n) = filter2(h, im(:,:,n));
    end
    
end