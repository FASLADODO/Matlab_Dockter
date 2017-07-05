% function [ edge_image ] = detect_edges( input_image , hsize , sigma , Thigh , Tlow)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


input_image = imread('bear.jpg');
hsize = 5;
sigma = 2;
MAX_Gdir = 1.5879;
Thigh = 0.3;%0.8*MAX_Gdir;
Tlow = 0.4*Thigh;

input_image = rgb2gray(input_image);

input_image = im2double(input_image);

% STEP 1: Convolve with the gaussian
gauss_kernel = fspecial('gaussian', hsize, sigma);
im_gauss = conv2(input_image,gauss_kernel,'same');

% STEP 2: calculate gradient intensity and direction
[Gmag,Gdir] = imgradient(im_gauss,'sobel');

figure;
imshow(Gmag);
title('Image Gradient');

% bin image direction 0, 45, 90, 135
Gdir = abs(Gdir);
Gdir(Gdir < 45) = 0;
Gdir(Gdir >= 45 & Gdir < 90) = 45;
Gdir(Gdir >= 90 & Gdir < 135) = 90;
Gdir(Gdir >= 135) = 135;


% STEP 3: Apply non-maximum suppression
nsize = 1; % kernel size each side
win_size = (2*nsize+1);
[m,n] = size(Gmag);
im_supp = zeros(m,n);
for i = nsize+1:m-nsize
    for j = nsize+1:n-nsize
        window = Gmag( i-nsize:i+nsize , j-nsize:j+nsize );
        angle = Gdir(i,j);
        
        % apply non-maximum suppression to everything not perpendicular to
        % angle
        switch(angle)
            case 0  % horizontal middle row
%                 window_vals = [window(2,1) window(2,2) window(2,3)];
                window_vals = window(nsize+1 , : );
            case 45 % off diagonal
%                 window_vals = [window(3,1) window(2,2) window(1,3)];
                window_vals = window( win_size:win_size-1:end-1);
            case 90 % vertical middle column
%                 window_vals = [window(1,2) window(2,2) window(3,2)];
                window_vals = window( : , nsize+1 );
            case 135 % diagonal
%                 window_vals = [window(1,1) window(2,2) window(3,3)];
                window_vals = window(1:win_size+1:end);

            otherwise
                disp('Invalid');
        end
        % get max value within the window
        max_window = max(max(window_vals));
        if (Gmag(i,j) == max_window)
            value = Gmag(i,j);
        else
            value = 0;
        end
        
        % copy the value at (i,j) to output if it is equal to the max
        im_supp(i,j) = value;
        
    end
end



figure;
imshow(im_supp);
title('Maxima-Suppresion');

% STEP 4: Apply Hysteresis
% two way thresholding\
im_hyst = im_supp;
im_hyst2 = im_supp;
im_hyst2(im_hyst2<Thigh) = 0;

im_hyst(im_hyst < Tlow) = 0; % suppress all values low values
for i = 2:m-1
    for j = 2:n-1
        % look for pixels between the high and low threshold
        if ( im_hyst(i,j) > Tlow && im_hyst(i,j) < Thigh )
            % see if any of these pixels touch a high edge
            % Check to see if there are any neighbors
            % using an 8 connected component
            ker_hyst = im_hyst(i-1:i+1 , j-1:j+1);
            ker_max = max(max(ker_hyst));
            if ( ker_max <= Thigh ) % if not touching high edge suppress         
                im_hyst(i,j) = 0;
            end
        end
    end
end

figure; imshow(im_hyst2);
title('Hysteresis2');

figure; imshow(im_hyst);
title('Hysteresis');

im_edge = edge(input_image);
figure; imshow(im_edge);
title('MATLAB Edge Function');
% end

