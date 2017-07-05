function bw = ThresholdImage(image,threshold,inv)
%image grayscale uint8
%threshold 0-255
% inv = 0 for makes image> threshold = 255
%inv = 1 makes image < threshold = 255;

if(inv)
    bwlog = image < threshold;
    bw = uint8(bwlog)*255;
else
    bwlog = image > threshold;
    bw = uint8(bwlog)*255;
end

end