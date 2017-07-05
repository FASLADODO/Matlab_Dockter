input = imread('eifel.jpg');

figure, imshow(input)

scale = [2,1.5];
zoomim = ZoomImageCenter(input,scale);


figure, imshow(zoomim)