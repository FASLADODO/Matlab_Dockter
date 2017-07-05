m = imread('Monkey.JPG');

g = rgb2gray(m);

imshow(m)

imshow(g)

d = im2double(g);

surf(d)

z = m(:,:,1);
